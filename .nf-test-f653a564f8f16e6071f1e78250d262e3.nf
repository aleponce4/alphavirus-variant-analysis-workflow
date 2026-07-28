import groovy.json.JsonGenerator
import groovy.json.JsonGenerator.Converter

nextflow.enable.dsl=2

// comes from nf-test to store json files
params.nf_test_output  = ""

// include dependencies


// include test process
include { IVAR_VARIANTS } from '/mnt/c/Users/aponcefl/alphavirus-variant-analysis-workflow/modules/local/ivar/variants/main.nf'

// define custom rules for JSON that will be generated.
def jsonOutput =
    new JsonGenerator.Options()
        .addConverter(Path) { value -> value.toAbsolutePath().toString() } // Custom converter for Path. Only filename
        .build()

def jsonWorkflowOutput = new JsonGenerator.Options().excludeNulls().build()


workflow {

    // run dependencies
    

    // process mapping
    def input = []
    
                input[0] = [
                    [ id: 'sampleA' ],
                    file("/mnt/c/Users/aponcefl/alphavirus-variant-analysis-workflow/tests/data/sampleA.test.bam", checkIfExists: true),
                    file("/mnt/c/Users/aponcefl/alphavirus-variant-analysis-workflow/tests/data/sampleA.test.bam.bai", checkIfExists: true)
                ]
                input[1] = [
                    [ id: 'viral_ref' ],
                    file("/mnt/c/Users/aponcefl/alphavirus-variant-analysis-workflow/tests/data/viral_ref.test.fasta", checkIfExists: true),
                    file("/mnt/c/Users/aponcefl/alphavirus-variant-analysis-workflow/tests/data/viral_ref.test.fasta.fai", checkIfExists: true)
                ]
                input[2] = [
                    [ id: 'viral_gff' ],
                    file("/mnt/c/Users/aponcefl/alphavirus-variant-analysis-workflow/tests/data/viral_ref.test.gff3", checkIfExists: true)
                ]
                
    //----

    //run process
    IVAR_VARIANTS(*input)

    if (IVAR_VARIANTS.output){

        // consumes all named output channels and stores items in a json file
        for (def name in IVAR_VARIANTS.out.getNames()) {
            serializeChannel(name, IVAR_VARIANTS.out.getProperty(name), jsonOutput)
        }	  
      
        // consumes all unnamed output channels and stores items in a json file
        def array = IVAR_VARIANTS.out as Object[]
        for (def i = 0; i < array.length ; i++) {
            serializeChannel(i, array[i], jsonOutput)
        }    	

    }
  
}

def serializeChannel(name, channel, jsonOutput) {
    def _name = name
    def list = [ ]
    channel.subscribe(
        onNext: {
            list.add(it)
        },
        onComplete: {
              def map = new HashMap()
              map[_name] = list
              def filename = "${params.nf_test_output}/output_${_name}.json"
              new File(filename).text = jsonOutput.toJson(map)		  		
        } 
    )
}


workflow.onComplete {

    def result = [
        success: workflow.success,
        exitStatus: workflow.exitStatus,
        errorMessage: workflow.errorMessage,
        errorReport: workflow.errorReport
    ]
    new File("${params.nf_test_output}/workflow.json").text = jsonWorkflowOutput.toJson(result)
    
}
