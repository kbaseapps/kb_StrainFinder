/*
A KBase module: kb_StrainFinder
*/

module kb_StrainFinder {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    typedef structure {
	string workspace_name;
        string in_genome_ref;
        /*string in_vcf_refs;*/  /* not ready yet */
        list<string> in_readslib_refs;
	string out_genomeSet_obj_name;

	int min_mapping_quality;
	int min_depth;
    } StrainFinder_v1_InputType;

    /*
        This example function accepts any number of parameters and returns results in a KBaseReport
    */
    funcdef run_StrainFinder_v1(StrainFinder_v1_InputType params) returns (ReportResults output) authentication required;

};
