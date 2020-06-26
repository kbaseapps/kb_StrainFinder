# -*- coding: utf-8 -*-
import os
import time
import unittest
import shutil
from pprint import pprint
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from kb_StrainFinder.kb_StrainFinderImpl import kb_StrainFinder
from kb_StrainFinder.kb_StrainFinderServer import MethodContext
from kb_StrainFinder.authclient import KBaseAuth as _KBaseAuth

from installed_clients.WorkspaceClient import Workspace
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.ReadsUtilsClient import ReadsUtils


class kb_StrainFinderTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_StrainFinder'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_StrainFinder',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = kb_StrainFinder(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_ContigFilter_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')
        if hasattr(cls, 'shock_ids'):
            for shock_id in cls.shock_ids:
                print('Deleting SHOCK node: '+str(shock_id))
                cls.delete_shock_node(shock_id)

    @classmethod
    def delete_shock_node(cls, node_id):
        header = {'Authorization': 'Oauth {0}'.format(cls.token)}
        requests.delete(cls.shockURL + '/node/' + node_id, headers=header,
                        allow_redirects=True)
        print('Deleted shock node ' + node_id)

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_kb_ReadsUtilities_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx


    # test_run_StrainFinder_v1_01
    # HIDE @unittest.skip("skipped test_run_StrainFinder_v1_01()")  # uncomment to skip

    def test_run_StrainFinder_v1_01(self):
        method = 'test_run_StrainFinder_v1_01'
        
        print ("\n\nRUNNING "+method+"()")
        print ("===========================================\n\n")

        # connect to clients
        try:
            gfuClient = GenomeFileUtil(self.callback_url, token=self.getContext()['token'])
        except Exception as e:
            raise ValueError('Unable to instantiate gfuClient with callbackURL: '+ self.callback_url +' ERROR: ' + str(e))
        try:
            ruClient = ReadsUtils(self.callback_url, token=self.getContext()['token'])
        except Exception as e:
            raise ValueError('Unable to instantiate ruClient with callbackURL: '+ self.callback_url +' ERROR: ' + str(e))

        # upload data
        sci_name = 'Thermodesulfobacterium thermophilum DSM 1276',
        base_genome = 'GCF_000421605.1_ASM42160v1_genomic'
        base_reads = 'seven_species_nonuniform_10K.PE_reads_paired-0'
        genome_gff_file = base_genome+'.gff.gz'
        genome_fna_file = base_genome+'.fna.gz'
        reads_file =  base_reads+'.fq'
        genome_gff_path = os.path.join(self.scratch, genome_gff_file)
        genome_fna_path = os.path.join(self.scratch, genome_fna_file)
        reads_path = os.path.join(self.scratch, reads_file)
        shutil.copy(os.path.join("data", genome_gff_file), genome_gff_path)
        shutil.copy(os.path.join("data", genome_fna_file), genome_fna_path)
        shutil.copy(os.path.join("data", reads_file), reads_path)
        genome_ref = gfuClient.fasta_gff_to_genome({
            'workspace_name': self.getWsName(),
            'fasta_file': {'path': genome_fna_path},
            'gff_file': {'path': genome_gff_path},
            'generate_missing_genes': 1,
            'source': 'GFF',
            #'scientific_name': sci_name,  # this is causing an error for some reason
            'genome_name': base_genome+'.Genome'
        }).get('genome_ref')

        reads_ref = ruClient.upload_reads({
            'wsname': self.getWsName(),
            'fwd_file': reads_path,
            'sequencing_tech': 'artificial reads',
            'interleaved': 1,
            'name': base_reads+'.Reads'
        })['obj_ref']
        
        # run method
        output_name = method+'_output'
        params = {
            'workspace_name': self.getWsName(),
            'in_genome_ref': genome_ref,
            'in_readslib_ref': reads_ref,
            'out_genomeSet_obj_name': output_name,
            'min_mapping_quality': 30,
            'min_depth': 50
        }
        result = self.getImpl().run_StrainFinder_v1(self.getContext(),params)

        print("\nRESULT:")
        pprint(result)
        pass
