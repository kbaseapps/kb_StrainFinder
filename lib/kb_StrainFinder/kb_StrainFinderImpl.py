# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import re
import subprocess
import sys
import shutil
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
#from installed_clients.SetAPIServiceClient import SetAPI_Service  # if you want to use the service wizard
from installed_clients.SetAPIClient import SetAPI  # if you want to run as SDK_LOCAL

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
    VERSION = "0.1.0"
    GIT_URL = "https://github.com/kbaseapps/kb_StrainFinder"
    GIT_COMMIT_HASH = "7189d399ad37ad6a735dd3138e50be5c24c1609a"

    #BEGIN_CLASS_HEADER

    # binaries
    STRAINFINDER_v1_installdir  = "/kb/module/strainfinder"
    STRAINFINDER_v1_bin         = os.path.join(STRAINFINDER_v1_installdir, "strainFinder.py")
    STRAINFINDER_v1_RUN_FIT_bin = os.path.join(STRAINFINDER_v1_installdir, "example", "run_fit.py")
    #VCFTOOLS_bin                = "/usr/local/bin/vcftools"
    
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

    def _translate_nuc_to_prot_seq(self, nuc_seq=None, genetic_code=None, keep_stop=False, truncate_to_stop=False):
        if not genetic_code:
            genetic_code = '11'
        if genetic_code != '11':
            raise ValueError('Method _translate_nuc_to_prot_seq() only knows genetic code 11')

        stop_char = '*'
        nuc_seq = nuc_seq.upper()
        prot_seq = ''

        genetic_code_table = dict()
        genetic_code_table['11'] = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':stop_char, 'TAG':stop_char,
            'TGC':'C', 'TGT':'C', 'TGA':stop_char, 'TGG':'W'
            }
        if genetic_code not in genetic_code_table:
            raise ValueError ("genetic code '"+str(genetic_code)+"' not configured in genetic_code_table")

        prot_seq = ''.join([genetic_code_table[genetic_code].get(nuc_seq[3*i:3*i+3],'X') for i in range(len(nuc_seq)//3)])

        if truncate_to_stop:
            new_prot_seq = ''
            for c in prot_seq:
                new_prot_seq += c
                if c == stop_char:
                    break
            prot_seq = new_prot_seq
                    
        if prot_seq.endswith(stop_char) and not keep_stop:
            prot_seq = prot_seq.rstrip(stop_char)
            
        return prot_seq

    def _reverse_complement (self, nuc_seq):
        nuc_seq = nuc_seq.upper()
        complement = { 'G': 'C',
                       'C': 'G',
                       'A': 'T',
                       'T': 'A',
                       'U': 'A'
                     }
        rev_nuc_seq = ''
        for c in nuc_seq[::-1]:
            rev_nuc_seq += complement[c]
        return rev_nuc_seq
    
    def read_fasta_file (self, fasta_file):
        fasta_buf = dict()
        headers = dict()
        id_order = []
        console = []
        
        with open (fasta_file, 'r') as fasta_handle:
            last_header = None
            last_id = None
            seq = ''
            for fasta_row in fasta_handle.readlines():
                fasta_row = fasta_row.rstrip()
                if fasta_row.startswith('>'):
                    this_header = fasta_row
                    this_id = this_header.replace('>','',1)
                    this_id = re.sub(' .*$','',this_id)
                    id_order.append(this_id)
                    if last_id != None:
                        headers[last_id] = last_header
                        fasta_buf[last_id] = seq
                    last_header = this_header
                    last_id = this_id
                    seq = ''
                else:
                    seq += fasta_row.replace(' ','')
            if last_id != None:
                headers[last_id] = last_header
                fasta_buf[last_id] = seq
                last_header = None
                last_id = None
                seq = ''

        return {'fasta': fasta_buf,
                'headers': headers,
                'id_order': id_order
                }
                
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

        self.SE_flag = 'SE'
        self.PE_flag = 'PE'

        #END_CONSTRUCTOR
        pass


    def run_StrainFinder_v1(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of type "StrainFinder_v1_InputType" ->
           structure: parameter "workspace_name" of String, parameter
           "in_genome_ref" of String, parameter "in_readslib_refs" of list of
           String, parameter "out_genomeSet_obj_name" of String, parameter
           "min_mapping_quality" of Long, parameter "min_depth" of Long
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_StrainFinder_v1

        #### STEP 0: Init
        ##
        DEBUG_MODE = 0
        method = 'run_StrainFinder_v1'
        console = []
        report_text = ''
        html_links = []
        file_links = []
        objects_created = []
        self.log(console, 'Running run_StrainFinder_v1() with parameters: ')
        self.log(console, "\n"+pformat(params))
        
        token = ctx['token']
        headers = {'Authorization': 'OAuth '+token}
        env = os.environ.copy()
        env['KB_AUTH_TOKEN'] = token
        
        output_dir = os.path.join(self.scratch, 'output_'+str(uuid.uuid4()))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # object_info tuple
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)

        # Client handles
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
        try:
            gfuClient = GenomeFileUtil (url=self.callbackURL, token=token)  # SDK local
        except Exception as e:
            raise ValueError('Unable to get GenomeFileUtil Client' +"\n" + str(e))
        try:
            auClient = AssemblyUtil (url=self.callbackURL, token=token)  # SDK local
        except Exception as e:
            raise ValueError('Unable to get AssemblyUtil Client' +"\n" + str(e))
        try:
            #setAPI_Client = SetAPI_Service (url=self.serviceWizardURL, token=token, service_ver=SERVICE_VER)  # Service
            setAPI_Client = SetAPI (url=self.callbackURL, token=token, service_ver=SERVICE_VER)  # SDK Local
        except Exception as e:
            raise ValueError('Unable to get SetAPI Client' +"\n" + str(e))
            

        #### STEP 1: Param checks
        ##
        required_params = ['workspace_name',
                           'in_genome_ref',
                           'in_readslib_refs',
                           'min_mapping_quality',
                           'min_depth',
                           'out_genomeSet_obj_name'
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

        
        #### STEP 2: Configure overall provenance
        ##
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        provenance[0]['input_ws_objects']=[]
        provenance[0]['input_ws_objects'].append(params['in_genome_ref'])
        provenance[0]['input_ws_objects'].extend(params['in_readslib_refs'])
        provenance[0]['service'] = 'kb_StrainFinder'
        provenance[0]['method'] = method
        

        #### STEP 3: Get Reads refs
        ##
        expanded_reads = []
        input_ref_seen = dict()
        SE_types = ['KBaseFile.SingleEndLibrary', 'KBaseAssembly.SingleEndLibrary']
        PE_types = ['KBaseFile.PairedEndLibrary', 'KBaseAssembly.PairedEndLibrary']

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple
        for input_ref in params['in_readslib_refs']:
            input_info = wsClient.get_object_info3({'objects': [{'ref': input_ref}]})['infos'][0]
            obj_name = input_info[NAME_I]
            type_name = input_info[TYPE_I].split('-')[0]

            # ReadsSet
            if type_name in ['KBaseSets.ReadsSet']:
                try:
                    input_readsSet_obj = setAPI_Client.get_reads_set_v1 ({'ref':input_ref,'include_item_info':1})
                except Exception as e:
                    raise ValueError('SetAPI FAILURE: Unable to get read library set object from workspace: (' + str(input_ref)+")\n" + str(e))

                for readsLibrary_obj in input_readsSet_obj['data']['items']:
                    this_reads_ref = readsLibrary_obj['ref']
                    if this_reads_ref in input_ref_seen:
                        continue
                    input_ref_seen[this_reads_ref] = True

                    this_reads_name = readsLibrary_obj['info'][NAME_I]
                    reads_item_type = readsLibrary_obj['info'][TYPE_I]
                    reads_item_type = re.sub ('-[0-9]+\.[0-9]+$', "", reads_item_type)  # remove trailing version
                    if reads_item_type in PE_types:
                        this_reads_type = self.PE_flag
                    elif reads_item_type in SE_types:
                        this_reads_type = self.SE_flag
                    else:
                        raise ValueError ("Can't handle read item type '"+reads_item_type+"' obj_name: '"+this_reads_name+" in Set: '"+str(input_ref)+"'")
                    expanded_reads.append({'ref':  this_reads_ref,
                                           'name': this_reads_name,
                                           'type': this_reads_type
                                          })
            # SingleEnd Library
            elif type_name in SE_types:
                this_reads_ref = input_ref
                if this_reads_ref in input_ref_seen:
                    continue
                input_ref_seen[this_reads_ref] = True
                this_reads_name = obj_name
                this_reads_type = self.SE_flag
                expanded_reads.append({'ref':  this_reads_ref,
                                       'name': this_reads_name,
                                       'type': this_reads_type
                                      })
            # PairedEnd Library
            elif type_name in PE_types:
                this_reads_ref = input_ref
                if this_reads_ref in input_ref_seen:
                    continue
                input_ref_seen[this_reads_ref] = True
                this_reads_name = obj_name
                this_reads_type = self.PE_flag
                expanded_reads.append({'ref':  this_reads_ref,
                                       'name': this_reads_name,
                                       'type': this_reads_type
                                      })
            else:
                raise ValueError ("Illegal type in input_refs: "+str(obj_name)+" ("+str(input_ref)+") is of type: '"+str(type_name)+"'")

        expanded_reads_refs = []
        for reads_element in expanded_reads:
            expanded_reads_refs.append(reads_element['ref'])
            
        """
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
        
        #### STEP 4: Get Genome's Assembly ref
        ##
        self.log(console, "GETTING GENOME ASSEMBLY OBJECT")
        try:
            genome_object_ret = wsClient.get_objects2({'objects':[{'ref':params['in_genome_ref']}]})['data'][0]
        except Exception as e:
            raise ValueError ("unable to get genome object "+params['in_genome_ref']+". "+str(e))
        genome_object = genome_object_ret['data']
        genome_obj_info = genome_object_ret['info']
        input_genome_obj_name = genome_obj_info[NAME_I]
        
        if not genome_object.get('assembly_ref'):
            raise ValueError ('OLD Genome type cannot be used with method.  Please update your Genome '+params['in_genome_ref']+' to modern Geomee object format')
        assembly_ref = genome_object['assembly_ref']
        
        
        #### STEP 5: Do Read mapping and get VCF variant info
        ##
        self.log(console, "READ MAPPING PHASE")
        ws_info = wsClient.get_workspace_info({'workspace': params['workspace_name']})
        workspace_id = ws_info[0]
        sub_method = 'call_variants'
        call_variants_params = {
            'workspace_name': params['workspace_name'],
            'workspace_id': workspace_id,
            'assembly_ref': assembly_ref,
            #'reads_ref': params['in_readslib_ref'],
            'reads_refs': expanded_reads_refs,
            'min_mapping_quality': params['min_mapping_quality'],
            'min_depth': params['min_depth'],
            'output_vcf': params['out_genomeSet_obj_name']+'.VCF'
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


        # save this report obj for later
        metadecoder_call_variants_reportObj = this_report_obj
        #self.log(console, pformat(metadecoder_call_variants_reportObj))  # DEBUG

        
        #### STEP 6: Extract VCFs
        ##
        ## Note: this should ultimately come from a VCF object with an API to spit it to file
        ##       currently will pull VCF from SHOCK directly.
        ##       from https://github.com/kbaseapps/kb_meta_decoder/blob/master/lib/kb_meta_decoder/kb_meta_decoderImpl.py
        ##       file_links[0] is BAM
        ##       file_links[1] is VCF
        ##
        """  It'll be something like this
        vcf_types = ['KBaseVariation.VCF']
        for obj in metadecoder_call_variants_reportObj['objects_created']: 
            obj_info = wsClient.get_obj_info_new({'objects':[{'ref':obj_ref}]})[0]
            obj_type = re.sub('-[0-9]+\.[0-9]+$', "", obj_info[TYPE_I])  # remove trailing version
            if obj_type in vcf_types:
                vcf_ref = obj_ref
                break
        """

        vcf_files = []
        for reads_lib_i,reads_ref in enumerate(expanded_reads_refs):
            self.log(console, "EXTRACTING VCF for ReadsLib "+str(reads_lib_i+1))
            vcf_shock_handle_id = None
            vcf_seen_cnt = 0;
            for file_link in metadecoder_call_variants_reportObj['file_links']:
                if re.search('VCF', file_link['label']):
                    vcf_seen_cnt += 1
                    if vcf_seen_cnt == (reads_lib_i+1):
                        vcf_shock_handle_id = file_link['handle']  # not 'handle_id'
                        break
            if not vcf_shock_handle_id:
                raise ValueError ("Failure to find VCF in SHOCK from alignment submethod")
            vcf_file = os.path.join(self.scratch, params['out_genomeSet_obj_name']+'.vcf')
            vcf_dl_result = dfuClient.shock_to_file({'handle_id': vcf_shock_handle_id,
                                                     'file_path': vcf_file,
                                                     'unpack': 'uncompress'
            })
            vcf_files.append(vcf_file)
            # DEBUG
            #with open(vcf_file, 'r') as vcf_file_handle:
            #    for line in vcf_file_handle.readlines():
            #        self.log(console, 'VCF_line: '+line)
            
        
        #### STEP 7: Parse VCF to get polymorphism frequencies
        ##
        run_dirs_list = []
        if DEBUG_MODE == 1:
            vcf_files = ["/kb/module/dev_test/data/Bin.002-37A_testACGT_ploidy1_alternate_2.vcf"]  # DEBUG]
        for reads_lib_i,reads_ref in enumerate(expanded_reads_refs):
            self.log(console, "Parsing VCF to polymorphism frequencies for ReadsLib "+str(reads_lib_i+1))
            vcf_file = vcf_files[reads_lib_i]
            vcf_buf = []
            SNP_freqs = dict()
                         
            with open (vcf_file, 'r') as vcf_handle:
                vcf_buf = vcf_handle.readlines()

            (CONTIG_ID_I, POS_I, SNP_ID_I, REF_SEQ_I, ALT_SEQ_I, QUAL_I, FILTER_I, INFO_I, FORMAT_I, EXTRA_INFO_I) = range(10)
            bases = ['A', 'C', 'G', 'T']
            base_i = { 'A': 0,
                       'C': 1,
                       'G': 2,
                       'T': 3,
                       'U': 3
            }
            for vcf_line in vcf_buf:
                if vcf_line.startswith('#'):
                    continue
                vcf_line = vcf_line.rstrip()
                row = vcf_line.split()
                contig_id      = row[CONTIG_ID_I]
                pos            = row[POS_I]
                ref_seq        = row[REF_SEQ_I]
                alt_seqs       = row[ALT_SEQ_I]
                var_info       = row[INFO_I]
                var_format     = row[FORMAT_I]
                var_extra_info = row[EXTRA_INFO_I]
            
                # just record SNPs
                if not var_format.startswith('GT:PL'):  # skip invariant positions
                    continue
                if var_info.startswith('INDEL'):  # can't handle INDELs yet
                    continue

                # init
                if contig_id not in SNP_freqs:
                    SNP_freqs[contig_id] = dict()
                if pos not in SNP_freqs[contig_id]:
                    SNP_freqs[contig_id][pos] = []
                    for i in range(4):
                        SNP_freqs[contig_id][pos].append(0)

                # get counts
                counts = []
                alt_seq_list = alt_seqs.split(',')
                [GT, PL, AD_counts_str] = var_extra_info.split(':')
                counts = list(map(int, AD_counts_str.split(',')))
                total_counts = 0
                for cnt in counts:
                    total_counts += cnt
                SNP_freqs[contig_id][pos][base_i[ref_seq]] = int(0.5 + 100.0 * counts[0] / float(total_counts))
                for alt_i,alt_seq in enumerate(alt_seq_list):
                    SNP_freqs[contig_id][pos][base_i[alt_seq]] = int(0.5 + 100.0 * float(counts[1+alt_i]) / float(total_counts))

                
            # write SNPs in strainfinder format and record position for later contig sequence
            #run_dir = os.path.join(self.STRAINFINDER_v1_installdir, 'example')
            this_run_dir = os.path.join(self.STRAINFINDER_v1_installdir, 'run_'+str(reads_lib_i))
            if not os.path.exists(this_run_dir):
                os.makedirs (this_run_dir)
            run_dirs_list.append(this_run_dir)
            allele_counts_file = os.path.join(this_run_dir, 'allele_counts.txt')
            if os.path.exists(allele_counts_file):
                shutil.move(allele_counts_file, allele_counts_file+'.orig-'+str(reads_lib_i))

            allele_counts_buf = []
            position_row_index = []
            allele_counts_buf.append("\t".join(['# A', 'C', 'G', 'T'])+"\n")
            for contig_id in sorted(SNP_freqs.keys()):
                for pos in sorted(SNP_freqs[contig_id].keys()):
                    position_row_index.append("\t".join([contig_id, pos]))
                    allele_counts_buf.append("\t".join(list(map(str, SNP_freqs[contig_id][pos])))+"\n")

            with open (allele_counts_file, 'w') as allele_counts_handle:
                allele_counts_handle.writelines(allele_counts_buf)

            # DEBUG
            #with open(allele_counts_file, 'r') as file_handle:
            #    for line in file_handle.readlines():
            #        self.log(console, 'ALLELE_line: '+line)

                
        #### STEP 8: Run StrainFinder
        ##
        fitted_genomes_files = []
        abund_vecs = []
        for reads_lib_i,reads_lib_ref in enumerate(expanded_reads_refs):
            self.log(console, "RUNNING STRAINFINDER for ReadsLib "+str(reads_lib_i+1))
            this_run_dir = run_dirs_list[reads_lib_i]
            fitted_genomes_file = os.path.join(this_run_dir, 'fitted_genomes.txt')
            if os.path.exists(fitted_genomes_file):
                os.remove(fitted_genomes_file)

            # Some subprocesses require shell=True in order to see input data
            #  also, if you do a redirect to out, you must join command first
            this_STRAINFINDER_v1_RUN_FIT_bin = os.path.join(this_run_dir, "run_fit.py")
            #shutil.copy(self.STRAINFINDER_v1_RUN_FIT_bin, this_STRAINFINDER_v1_RUN_FIT_bin)
            bin_buf = []
            with open (self.STRAINFINDER_v1_RUN_FIT_bin, 'r') as src_bin_handle:
                for line in src_bin_handle.readlines():
                    if 'True strain relative abundances' in line:
                        continue
                    bin_buf.append(line)
            with open (this_STRAINFINDER_v1_RUN_FIT_bin, 'w') as dst_bin_handle:
                dst_bin_handle.writelines(bin_buf)

            strainfinder_cmd = ['python']
            strainfinder_cmd.append(this_STRAINFINDER_v1_RUN_FIT_bin)
        
            #env = os.environ.copy()
            #p = subprocess.Popen([joined_cmd], \
            p = subprocess.Popen(strainfinder_cmd, \
                                 cwd = this_run_dir, \
                                 stdout = subprocess.PIPE, \
                                 stderr = subprocess.STDOUT, \
                                 shell = False)
                                 #env = env)
        
            # Read output
            #
            while True:
                line = p.stdout.readline()
                #line = p.stderr.readline()
                if not line: break
                self.log(console, line.replace('\n', ''))
                if line.startswith('Inferred strain relative abundances = ['):
                    inferred_abundances = line.replace('Inferred strain relative abundances = [','')
                    inferred_abundances = inferred_abundances.replace(']','')
                    inferred_abundances = inferred_abundances.strip()
                    inferred_abundances = inferred_abundances.replace('  ',' ')
                    abund_vec = inferred_abundances.split()
                    abund_vecs.append(abund_vec)
                    
            p.stdout.close()
            #p.stderr.close()
            p.wait()
            self.log(console, 'return code: ' + str(p.returncode))
            if p.returncode != 0:
                raise ValueError('Error running STRAINFINDER, return code: '+str(p.returncode) + 
                                 '\n\n'+ '\n'.join(console))

            # Check that STRAINFINDER produced output
            #
            if not os.path.isfile(fitted_genomes_file):
                raise ValueError("failed to create STRAINFINDER output: "+fitted_genomes_file)
            fitted_genomes_files.append(fitted_genomes_file)
            

        #### STEP 9: download in genome assembly fasta file, gff file, and read fields from Genome obj
        ##
        ## Note: running call_variants() appears to delete the input assembly
        ##  so must download that assembly AFTER (or copy it to somewhere safe)
        ##  we're doing it now.
        src_scientific_name = None
        src_source = None
        src_release = None
        src_genetic_code = None
        #src_taxon_wsname = None
        #src_taxon_id = None

        # Get genome obj to read fields to set on GFF upload
        try:
            input_genome_obj = wsClient.get_objects2({'objects':[{'ref': params['in_genome_ref']}]})['data'][0]['data']
        except:
            raise ValueError ("unable to get genome obj "+params['in_genome_ref']+". "+str(e))

        if 'scientific_name' in input_genome_obj:
            src_scientific_name = input_genome_obj['scientific_name']
        if 'source' in input_genome_obj:
            src_source = input_genome_obj['source']
        if 'release' in input_genome_obj:
            src_release = input_genome_obj['release']
        if 'genetic_code' in input_genome_obj:
            src_genetic_code = input_genome_obj['genetic_code']
        
        # Get genome assembly as fasta file
        try:
            input_genome_fasta_file = auClient.get_assembly_as_fasta({'ref':assembly_ref})['path']
        except Exception as e:
            raise ValueError ("unable to get genome fasta "+params['in_genome_ref']+". "+str(e))
        # DEBUG
        if DEBUG_MODE == 1:
            input_genome_fasta_file = "/kb/module/dev_test/data/Bin.002.Genome.assembly.fa"

        fasta_read = self.read_fasta_file(input_genome_fasta_file)
        base_genome_fasta = fasta_read['fasta']
        base_genome_headers = fasta_read['headers']
        base_genome_contigID_order = fasta_read['id_order']

        # Get GFF as file
        try:
            input_genome_gff_file = gfuClient.genome_to_gff({
                'genome_ref': params['in_genome_ref'],
                'target_dir': output_dir})['file_path']
        except:
            raise ValueError ("unable to get genome gff "+params['in_genome_ref']+". "+str(e))
        # DEBUG
        if DEBUG_MODE == 1:
            input_genome_gff_file = "/kb/module/dev_test/data/Bin.002.fasta_assembly.RAST.gff"

            
        #### STEP 10: Create Assembly FASTAs for strain genomes
        ##
        #num_genomes_found_per_readslib = []
        all_new_genome_refs = []
        all_new_genome_names = []
        all_set_elements = dict()
        num_strain_genomes_generated = 0
        for reads_lib_i,reads_lib_ref in enumerate(expanded_reads_refs):
            self.log(console, "GETTING ALLELES FOR STRAIN MODES for ReadsLib "+str(reads_lib_i+1))
            fitted_genomes_file = fitted_genomes_files[reads_lib_i]
            fitted_genomes_rows = []
            with open (fitted_genomes_file, 'r') as fitted_genomes_handle:
                for row in fitted_genomes_handle.readlines():
                    fitted_genomes_rows.append(row.rstrip().split())
            num_genomes_found = len(fitted_genomes_rows[0])
            #num_genomes_found_per_readslib.append(num_genomes_found)
            
            # only one strain mode (may not be same as reference so do generate strain genome
            if num_genomes_found == 1:
                msg = "ReadsLib "+str(reads_lib_i)+": StrainFinder found only one strain in the data.  Not creating GenomeSet for Reads Library "+str(reads_lib_i+1)
                self.log(console, msg)
                report_text += msg+"\n"

            # multiple strain modes
            num_strain_genomes_generated += num_genomes_found
                
            # set provenance to just include this reads lib for input_ws_objects
            this_provenance = provenance
            this_provenance[0]['input_ws_objects']=[]
            this_provenance[0]['input_ws_objects'].append(params['in_genome_ref'])
            this_provenance[0]['input_ws_objects'].append(reads_lib_ref)
            
            ####  Make revised Assembly Fasta
            ##
            new_fasta_files = []
            for genome_i in range(num_genomes_found):
                new_genome_fasta = dict()
                new_genome_headers = dict()
                new_genome_fasta_len = dict()
                for contigID in base_genome_contigID_order:
                    new_genome_fasta[contigID] = base_genome_fasta[contigID]
                    new_genome_headers[contigID] = base_genome_headers[contigID]
                    new_genome_fasta_len[contigID] = len(base_genome_fasta[contigID])
                for row_i,row in enumerate(fitted_genomes_rows):
                    this_allele = row[genome_i]
                    [this_contigID, this_pos_n] = position_row_index[row_i].split("\t")
                    pos_i = int(this_pos_n)-1
                    if pos_i == 0:
                        new_genome_fasta[this_contigID] = this_allele + new_genome_fasta[this_contigID][pos_i+1:]
                    elif pos_i == new_genome_fasta_len[this_contigID]-1:
                        new_genome_fasta[this_contigID] = new_genome_fasta[this_contigID][:pos_i] + this_allele
                    else:
                        new_genome_fasta[this_contigID] = new_genome_fasta[this_contigID][:pos_i] + this_allele + new_genome_fasta[this_contigID][pos_i+1:]

                this_fasta_file = os.path.join(output_dir, 'new_fasta-'+str(reads_lib_i)+'-'+str(genome_i)+'-'+str(uuid.uuid4())+'.fasta')
                new_fasta_files.append(this_fasta_file)
                with open (this_fasta_file, 'w') as this_fasta_handle:
                    for contigID in base_genome_contigID_order:
                        this_fasta_handle.write(new_genome_headers[contigID]+"\n")
                        this_fasta_handle.write(new_genome_fasta[contigID]+"\n")
                
                        
            #### Check for stop codons and if yes, adjust GFFs
            ##
            self.log(console, "READING FEATURES")
            SNP_pos = dict()
            SNP_CDS = dict()

            # initialize and determine SNP positions
            for contigID in base_genome_contigID_order:
                SNP_pos[contigID] = []
                SNP_CDS[contigID] = dict()
            for row_i,row in enumerate(fitted_genomes_rows):
                [this_contigID, this_pos_n] = position_row_index[row_i].split("\t")
                SNP_pos[this_contigID].append(int(this_pos_n))

            # read base genome feature locations and determine which have SNPs
            with open (input_genome_gff_file, 'r') as gff_handle:
                for gff_line in gff_handle.readlines():
                    gff_line = gff_line.rstrip()
                    if gff_line.startswith('#'):
                        continue
                    [contigID, annot_db, locus_type, beg_str, end_str, d1, strand, d2, info] = gff_line.split("\t")
                    beg = int(beg_str)
                    end = int(end_str)
                    feature_has_SNP = False
                    for pos_n in SNP_pos[contigID]:
                        if pos_n >= beg and pos_n <= end:
                            feature_has_SNP = True
                            break
                    if feature_has_SNP:
                        report_text += gff_line+"\n"
                        if locus_type == 'CDS':
                            loc = ",".join([str(beg),str(end),strand])
                            SNP_CDS[contigID][loc] = True

            # build strain genomes
            new_gff_files = []
            for genome_i in range(num_genomes_found):
                new_gff_buf = []
                with open (input_genome_gff_file, 'r') as gff_handle:
                    for gff_line in gff_handle.readlines():
                        gff_line = gff_line.rstrip()
                        if gff_line.startswith('#'):
                            new_gff_buf.append(gff_line)
                            continue
                        [contigID, annot_db, locus_type, beg_str, end_str, d1, strand, d2, info] = gff_line.split("\t")
                        # this will cover both CDS and parent gene
                        loc = ",".join([str(beg),str(end),strand])
                        if loc not in SNP_CDS[contigID]:  
                            new_gff_buf.append(gff_line)
                            continue
                        # else check for altered stop
                        contig_len = len(base_genome_fasta[contigID])
                        if end <= contig_len:
                            nuc_seq = base_genome_fasta[contigID][beg-1:end]
                        else:
                            nuc_seq = base_genome_fasta[contigID][beg-1:contig_len] + \
                                      base_genome_fasta[contigID][0:end-contig_len+1]
                        if strand == '-':
                            nuc_seq = self._reverse_complement(nuc_seq)
                        prot_seq = self._translate_nuc_to_prot_seq(nuc_seq=nuc_seq, keep_stop=True, truncate_to_stop=True)
                        old_len = end - beg + 1
                        new_len = 3 * len(prot_seq)
                        if new_len == old_len and prot_seq.endswith('*'):
                            new_gff_buf.append(gff_line)
                            continue

                        if not prot_seq.endswith('*'):
                            if strand == '-':
                                new_nuc_seq = self._reverse_complement(base_genome_fasta[contigID][0:end])
                            else:
                                new_nuc_seq = base_genome_fasta[contigID][beg-1:contig_len]
                            new_prot_seq = self._translate_nuc_to_prot_seq(nuc_seq=new_nuc_seq, keep_stop=True, truncate_to_stop=True)
                            new_len = 3 * len(new_prot_seq)

                        # adjust feature coords
                        if strand == '-':
                            new_beg = end - new_len+1
                            new_end = end
                        else:
                            new_beg = beg 
                            new_end = beg + new_len-1

                        if new_end <= contig_len:
                            new_gff_line = "\t".join([contigID, annot_db, locus_type, str(new_beg), str(new_end), d1, strand, d2, info])
                            new_gff_buf.append(new_gff_line)
                        else:
                            if new_beg > contig_len:
                                new_beg_p1 = new_beg - contig_len
                                new_end_p1 = new_end - contig_len
                                new_gff_line = "\t".join([contigID, annot_db, locus_type, str(new_beg_p1), str(new_end_p1), d1, strand, d2, info])
                                new_gff_buf.append(new_gff_line)
                            else:  # break into two
                                new_beg_p1 = new_beg
                                new_end_p1 = contig_len
                                new_beg_p2 = 1
                                new_end_p2 = new_end - contig_len
                                new_gff_line_p1 = "\t".join([contigID, annot_db, locus_type, str(new_beg_p1), str(new_end_p1), d1, strand, d2, info])
                                new_gff_line_p2 = "\t".join([contigID, annot_db, locus_type, str(new_beg_p2), str(new_end_p2), d1, strand, d2, info])
                                new_gff_buf.append(new_gff_line_p1)
                                new_gff_buf.append(new_gff_line_p2)
                                    
                this_gff_file = os.path.join(output_dir, 'new_gff-'+str(reads_lib_i)+'-'+str(genome_i)+'-'+str(uuid.uuid4())+'.gff')
                new_gff_files.append(this_gff_file)
                with open (this_gff_file, 'w') as this_gff_handle:
                    for gff_line in new_gff_buf:
                        this_gff_handle.write(gff_line+"\n")

                        
            #### STEP 12: Upload Strain Genomes and make GenomeSet for this ReadsLib
            ##
            self.log(console, "UPLOADING GENOMES for ReadsLib "+str(reads_lib_i+1))
            new_genome_refs = []
            new_genome_names = []
            set_elements = dict()
            #items = []

            for genome_i in range(num_genomes_found):
                new_genome_obj_name = re.sub(r'\.[^\.]+$','',params['out_genomeSet_obj_name'])
                new_genome_obj_name += '-Reads_'+str(reads_lib_i+1)+'-Strain_'+str(genome_i+1)
                new_genome_obj_name += ".Genome"
                self.log(console, "UPLOADING "+new_genome_obj_name)
                new_genome_ref = gfuClient.fasta_gff_to_genome({
                    'workspace_name': params['workspace_name'],
                    'fasta_file': {'path': new_fasta_files[genome_i]},
                    'gff_file': {'path': new_gff_files[genome_i]},
                    'genome_name': new_genome_obj_name,
                    'scientific_name': src_scientific_name,
                    'source': src_source,
                    'release': src_release,
                    'genetic_code': src_genetic_code
                })['genome_ref']

                new_genome_names.append(new_genome_obj_name)
                new_genome_refs.append(new_genome_ref)
                all_new_genome_names.append(new_genome_obj_name)
                all_new_genome_refs.append(new_genome_ref)
                    
                set_elements[new_genome_obj_name] = dict()
                set_elements[new_genome_obj_name]['ref'] = new_genome_ref
                all_set_elements[new_genome_obj_name] = dict()
                all_set_elements[new_genome_obj_name]['ref'] = new_genome_ref

            # attach created Genome objs to report
            for genome_i in range(num_genomes_found):
                objects_created.append({'ref': new_genome_refs[genome_i],
                                            'description': new_genome_names[genome_i]+' StrainFinder Genome'})

            # DEBUG
            #self.log(console, 'SET:'+"\n"+pformat(set_elements))
            #self.log(console, 'PROV:'+"\n"+pformat(this_provenance))
            # HERE
                    
            # create GenomeSet for this ReadsLib
            if num_genomes_found > 1 and len(expanded_reads_refs) > 1:
                genomeSet_name = params['out_genomeSet_obj_name']+'-Reads_'+str(reads_lib_i+1)
                genomeSet_obj = { 'description': 'Strain Genomes of '+input_genome_obj_name+' from ReadsLib '+str(reads_lib_i+1),
                                  #'items': items
                                  'elements': set_elements
                }
                try:
                    #output_genomeSet_ref = setAPI_Client.save_genome_set_v1 ({'workspace_name': params['workspace_name'],
                    #                                                          'output_object_name': genomeSet_name,
                    #                                                          'data': genomeSet_obj
                    #                                                      })['set_ref']
                    new_obj_info = wsClient.save_objects({'workspace': params['workspace_name'],
                                                          'objects': [{'type': 'KBaseSearch.GenomeSet',
                                                                       'data': genomeSet_obj,
                                                                       'name': genomeSet_name,
                                                                       'meta': {},
                                                                       'provenance': this_provenance
                                                          }]
                    })[0]
                except Exception as e:
                    raise ValueError('SetAPI FAILURE: Unable to save genome set object to workspace: (' + params['workspace_name']+")\n" + str(e))
                genomeSet_ref = '/'.join([str(new_obj_info[WORKSPACE_I]),
                                          str(new_obj_info[OBJID_I]),
                                          str(new_obj_info[VERSION_I])])
                
                objects_created.append({'ref': genomeSet_ref,
                                        'description': 'StrainFinder GenomeSet for ReadsLib '+str(reads_lib_i+1)})


            #### STEP 13: Add haplotype impact on genotype to Report?  Just minimally for now
            ##
            self.log(console, "ANALYZING GENOTYPES")
            for genome_i in range(num_genomes_found):
                report_text += str(100*(float(abund_vec[genome_i])))+' %'+"\n"
                
                    
        #### STEP 14: create GenomeSet for all generated Strains
        ##
        if num_strain_genomes_generated > 1:
            genomeSet_name = params['out_genomeSet_obj_name']
            genomeSet_obj = { 'description': 'Strain Genomes of '+input_genome_obj_name+' from ALL ReadsLibs',
                              #'items': items
                              'elements': all_set_elements
            }
            try:
                #output_genomeSet_ref = setAPI_Client.save_genome_set_v1 ({'workspace_name': params['workspace_name'],
                #                                                          'output_object_name': genomeSet_name,
                #                                                          'data': genomeSet_obj
                #                                                      })['set_ref']
                new_obj_info = wsClient.save_objects({'workspace': params['workspace_name'],
                                                      'objects': [{'type': 'KBaseSearch.GenomeSet',
                                                                   'data': genomeSet_obj,
                                                                   'name': genomeSet_name,
                                                                   'meta': {},
                                                                   'provenance': provenance
                                                      }]
                })[0]
            except Exception as e:
                raise ValueError('SetAPI FAILURE: Unable to save genome set object to workspace: (' + params['workspace_name']+")\n" + str(e))
            genomeSet_ref = '/'.join([str(new_obj_info[WORKSPACE_I]),
                                      str(new_obj_info[OBJID_I]),
                                      str(new_obj_info[VERSION_I])])
            
            objects_created.append({'ref': genomeSet_ref,
                                    'description': 'StrainFinder GenomeSet for ALL ReadsLibs'})


        #### STEP 15: Create Report
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

        # attach created objects and files to final report
        #objects_created.append(metadecoder_call_variants_reportObj['objects_created'])
        for file_link_item in metadecoder_call_variants_reportObj['file_links']:
            #this_shock_id = file_link_item['URL']
            this_shock_id = re.sub('^.*/', '', file_link_item['URL'])
            new_file_link_item = {'shock_id': this_shock_id,
                                  'name': file_link_item['name'],
                                  'label': file_link_item['label']
            }
            file_links.append(new_file_link_item)

        for html_link_item in metadecoder_call_variants_reportObj['html_links']:
            #this_shock_id = html_link_item['URL']
            this_shock_id = re.sub('^.*/', '', html_link_item['URL'])
            new_html_link_item = {'shock_id': this_shock_id,
                                  'name': html_link_item['name'],
                                  'label': html_link_item['label']
            }
            html_links.append(new_html_link_item)

        # until we replace with StrainFinder report info
        reportObj['direct_html_link_index'] = metadecoder_call_variants_reportObj['direct_html_link_index']
        reportObj['html_links'] = html_links
        reportObj['file_links'] = file_links
        reportObj['objects_created'] = objects_created
        if report_text:
            reportObj['message'] = report_text
            

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
