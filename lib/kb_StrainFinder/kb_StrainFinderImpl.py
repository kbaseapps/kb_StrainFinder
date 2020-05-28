# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import re
import subprocess
import sys
import traceback
import uuid
from datetime import datetime
from pprint import pformat

from installed_clients.WorkspaceClient import Workspace as workspaceService

from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.ReadsUtilsClient import ReadsUtils
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.SetAPIClient import SetAPI

from installed_clients.kb_meta_decoderClient import kb_meta_decoder

#END_HEADER


class kb_StrainFinder:
    '''
    Module Name:
    kb_StrainFinder

    Module Description:
    A KBase module: kb_StrainFinder
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/kbaseapps/kb_StrainFinder"
    GIT_COMMIT_HASH = "a856ee03f3469c091570fd7742862a2b93ea0570"

    #BEGIN_CLASS_HEADER

    # binaries
    STRAINFINDER_v1_bin = "/kb/module/strainfinder/strainFinder.py"
    VCFTOOLS_bin        = "/usr/local/bin/vcftools"
    
    # timestamp
    def now_ISO(self):
        now_timestamp = datetime.now()
        now_secs_from_epoch = (now_timestamp - datetime(1970,1,1)).total_seconds()
        now_timestamp_in_iso = datetime.fromtimestamp(int(now_secs_from_epoch)).strftime('%Y-%m-%d_%T')
        return now_timestamp_in_iso

    # message logging
    def log(self, target, message):
        message = '['+self.now_ISO()+'] '+message
        if target is not None:
            target.append(message)
        print(message)
        sys.stdout.flush()

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.shared_folder = config['scratch']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        self.config = config
        self.workspaceURL = config['workspace-url']
        self.shockURL = config['shock-url']
        self.handleURL = config['handle-service-url']
        self.serviceWizardURL = config['srv-wiz-url']
        self.callbackURL = os.environ.get('SDK_CALLBACK_URL')
        if self.callbackURL == None:
            raise ValueError ("SDK_CALLBACK_URL not set in environment")

        self.scratch = os.path.abspath(config['scratch'])
        if self.scratch == None:
            self.scratch = os.path.join('/kb','module','local_scratch')
        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)
        #END_CONSTRUCTOR
        pass


    def run_StrainFinder_v1(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of type "StrainFinder_v1_InputType" ->
           structure: parameter "workspace_name" of String, parameter
           "in_genome_ref" of String, parameter "in_readslib_ref" of String,
           parameter "out_obj_name" of String
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_StrainFinder_v1

        #### Init
        ##
        method = 'run_StrainFinder_v1'
        console = []
        report = ''
        file_links = []
        objects_created = []
        self.log(console, 'Running KButil_Split_Reads() with parameters: ')
        self.log(console, "\n"+pformat(params))
        
        token = ctx['token']
        headers = {'Authorization': 'OAuth '+token}
        env = os.environ.copy()
        env['KB_AUTH_TOKEN'] = token
        
        # object_info tuple
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)

        # Client handles
        #
        #SERVICE_VER = 'dev'  # DEBUG
        SERVICE_VER = 'release'
        try:
            wsClient = workspaceService(self.workspaceURL, token=token)
        except Exception as e:
            raise ValueError('Unable to get Workspace Client' +"\n" + str(e))
        try:
            dfuClient = DataFileUtil (url=self.callbackURL, token=token)  # SDK local
        except Exception as e:
            raise ValueError('Unable to get DataFileUtil Client' +"\n" + str(e))
            

        #### Param checks
        ##
        required_params = ['workspace_name',
                           'in_genome_ref',
                           'in_readslib_ref',
                           'out_obj_name'
                           ]
        for required_param in required_params:
            if required_param not in params or params[required_param] == None:
                raise ValueError ("Must define required param: '"+required_param+"'")
            
        """
        # and param defaults
        defaults = { 'split_num': 10
                   }
        for arg in defaults.keys():
            if arg not in params or params[arg] == None or params[arg] == '':
                params[arg] = defaults[arg]
        """

        
        #### Configure provenance
        ##
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        provenance[0]['input_ws_objects']=[]
        provenance[0]['input_ws_objects'].append(params['in_genome_ref'])
        provenance[0]['input_ws_objects'].append(params['in_readslib_ref'])
        provenance[0]['service'] = 'kb_StrainFinder'
        provenance[0]['method'] = method
        

        """
        #### Get Reads
        ##
        PairedEndTypes = ["KBaseFile.PairedEndLibrary","KBaseAssembly.PairedEndLibrary"]
        SingleEndTypes = ["KBaseFile.SingleEndLibrary","KBaseAssembly.SingleEndLibrary"]
        acceptable_types = PairedEndTypes + SingleEndTypes
        try:            
            input_reads_ref = params['in_readslib_ref']
            input_reads_obj_info = wsClient.get_object_info_new ({'objects':[{'ref':input_reads_ref}]})[0]
            input_reads_obj_type = input_reads_obj_info[TYPE_I]
            input_reads_obj_type = re.sub ('-[0-9]+\.[0-9]+$', "", input_reads_obj_type)  # remove trailing version
            #input_reads_obj_version = input_reads_obj_info[VERSION_I]  # this is object version, not type version
        except Exception as e:
            raise ValueError('Unable to get read library object info from workspace: (' + 
str(input_reads_ref) +')' + str(e))
        if input_reads_obj_type not in acceptable_types:
            raise ValueError ("Input reads of type: '"+input_reads_obj_type+"'.  Must be one of "+", ".join(acceptable_types))        

        
        # Download Reads
        self.log (console, "DOWNLOADING READS")  # DEBUG
        try:
            readsUtils_Client = ReadsUtils (url=self.callbackURL, token=ctx['token'])  # SDK local
        except Exception as e:
            raise ValueError('Unable to get ReadsUtils Client' +"\n" + str(e))
        try:
            readsLibrary = readsUtils_Client.download_reads ({'read_libraries': [input_reads_ref],
                                                             'interleaved': 'true'
                                                             })
        except Exception as e:
            raise ValueError('Unable to download read library sequences from workspace: (' + str(input_reads_ref) +")\n" + str(e))

        if input_reads_obj_type in PairedEndTypes:
            # Download reads Libs to FASTQ files
            in_reads_file_path = readsLibrary['files'][input_reads_ref]['files']['fwd']
        else:  # if input_reads_obj_type in SingleEndTypes:
            in_reads_file_path = readsLibrary['files'][input_reads_ref]['files']['fwd']
            
        """

        
        #### Get Genome's Assembly ref
        ##
        self.log(console, "GETTING GENOME ASSEMBLY OBJECT")
        genome_object = wsClient.get_objects2({'objects':[{'ref':params['in_genome_ref']}]})['data'][0]['data']
        if not genome_object.get('assembly_ref'):
            raise ValueError ('OLD Genome type cannot be used with method.  Please update your Genome '+params['in_genome_ref']+' to modern Geomee object format')
        assembly_ref = genome_object['assembly_ref']


        #### Do Read mapping and get VCF variant info
        ##
        self.log(console, "READ MAPPING PHASE")
        ws_info = wsClient.get_workspace_info({'workspace': params['workspace_name']})
        workspace_id = ws_info[0]
        sub_method = 'call_variants'
        call_variants_params = {
            'workspace_name': params['workspace_name'],
            'workspace_id': workspace_id,
            'assembly_ref': assembly_ref,
            'reads_ref': params['in_readslib_ref'],
            'output_vcf': params['out_obj_name']+'.VCF'
        }
        #MD_SERVICE_VER = 'release'
        MD_SERVICE_VER = 'dev'
        try:
            mdClient = kb_meta_decoder(self.callbackURL, token=token, service_ver=MD_SERVICE_VER)
        except Exception as e:
            raise ValueError("unable to instantiate metadecoderClient. "+str(e))
        try:
            self.log(console, "RUNNING call_variants()")
            this_retVal = mdClient.call_variants(call_variants_params)
        except Exception as e:
            raise ValueError ("unable to run "+sub_method+". "+str(e))
        try:
            this_report_obj = wsClient.get_objects2({'objects':[{'ref':this_retVal['report_ref']}]})['data'][0]['data']
        except Exception as e:
            raise ValueError("unable to fetch "+sub_method+" report: " + this_retVal['report_ref']+". "+str(e))
        objects_created.extend(this_report_obj['objects_created'])

        # attach created objects and files to final report
        objects_created.extend(this_report_obj['objects_created'])
        for file_link_item in this_report_obj['file_links']:
            #this_shock_id = file_link_item['URL']
            this_shock_id = re.sub('^.*/', '', file_link_item['URL'])
            new_file_link_item = {'shock_id': this_shock_id,
                                  'name': file_link_item['name'],
                                  'label': file_link_item['label']
            }
            file_links.append(new_file_link_item)

        # save this report obj for later
        metadecoder_call_variants_reportObj = this_report_obj
        #self.log(console, pformat(metadecoder_call_variants_reportObj))  # DEBUG
        

        #### Extract VCF
        ##
        ## Note: this should ultimately come from a VCF object with an API to spit it to file
        ##       currently will pull VCF from SHOCK directly.
        ##       from https://github.com/kbaseapps/kb_meta_decoder/blob/master/lib/kb_meta_decoder/kb_meta_decoderImpl.py
        ##       file_links[0] is BAM
        ##       file_links[1] is VCF
        ##
        self.log(console, "EXTRACTING VCF")
        """  It'll be something like this
        vcf_types = ['KBaseVariation.VCF']
        for obj in metadecoder_call_variants_reportObj['objects_created']: 
            obj_info = wsClient.get_obj_info_new({'objects':[{'ref':obj_ref}]})[0]
            obj_type = re.sub('-[0-9]+\.[0-9]+$', "", obj_info[TYPE_I])  # remove trailing version
            if obj_type in vcf_types:
                vcf_ref = obj_ref
                break
        """
        vcf_shock_handle_id = None
        for file_link in metadecoder_call_variants_reportObj['file_links']:
            if re.search('VCF', file_link['label']):
                vcf_shock_handle_id = file_link['handle']  # not 'handle_id'
                break
        if not vcf_shock_handle_id:
            raise ValueError ("Failure to find VCF in SHOCK from alignment submethod")
        vcf_file = os.path.join(self.scratch, params['out_obj_name']+'.vcf')
        vcf_dl_result = dfuClient.shock_to_file({'handle_id': vcf_shock_handle_id,
                                                 'file_path': vcf_file,
                                                 'unpack': 'uncompress'
                                                 })
        # DEBUG
        #with open(vcf_file, 'r') as vcf_file_handle:
        #    for line in vcf_file_handle.readlines():
        #        self.log(console, 'VCF_line: '+line)
                
        
        #### Create Report
        ##
        self.log(console, "CREATING REPORT")

        # instantiate report object
        reportName = method+'_report_' + str(uuid.uuid4())
        reportObj = {'objects_created': [],
                     'direct_html_link_index': 0,
                     'file_links': [],
                     'html_links': [],
                     'workspace_name': params['workspace_name'],
                     'report_object_name': reportName
                     }
        # can't just copy substructures because format of those fields in report
        #  object different from the format needed to pass to create_extended_report() method.
        #  for example, below doesn't work
        #for field in ('direct_html_link_index', 'file_links', 'html_links'):
        #    reportObj[field] = view_tree_reportObj[field]
        #    self.log<(console, "REPORT "+field+": "+pformat(view_tree_reportObj[field]))  # DEBUG
        #
        reportObj['direct_html_link_index'] = metadecoder_call_variants_reportObj['direct_html_link_index']
        for html_link_item in metadecoder_call_variants_reportObj['html_links']:
            #this_shock_id = html_link_item['URL']
            this_shock_id = re.sub('^.*/', '', html_link_item['URL'])
            new_html_link_item = {'shock_id': this_shock_id,
                                  'name': html_link_item['name'],
                                  'label': html_link_item['label']
            }
            reportObj['html_links'].append(new_html_link_item)

        reportObj['file_links'] = file_links
        reportObj['objects_created'] = objects_created


        # save report object
        try:
            reportClient = KBaseReport(self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
        except:
            raise ValueError ("unable to instantiate KBaseReport")
        report_info = reportClient.create_extended_report(reportObj)

        # Done
        self.log(console, "BUILDING RETURN OBJECT")
        output = {'report_name': report_info['name'],
                  'report_ref': report_info['ref']
        }
        #END run_StrainFinder_v1

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_StrainFinder_v1 return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
