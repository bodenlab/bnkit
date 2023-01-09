# JSON API 

## Requests for `asr.GServer`

    { "Command":<request>,
      "Job":<job-number> }

### `<request>` is `"Retrieve"`
Recall the request

### `<request>` is `"Output"`
Return the result of the request

    { "Job":<job-number>,
      "Result":{<result-JSON>} }

### `<request>` is `"Place"`
Find the place in queue of the request


## Commands for `asr.GServer`
Execute/submit job

    { "Command":<command>,
      "Auth":<auth-token>,
      "Params":<params> }

Server returns (if queued) 
    
    { "Message":"Queued","Job":<job-number> }

or runs job directly (based on `request.isQueued()`) and returns result

### Create input for reconstruction: `<command>` is `"Pogit"`

`<params>` is

    { "Tree":<tree>,
      "Alignment":<alignment> }

Result is a `<POGTree>`

    { 
      }

### Run reconstruction: `<command>` is `"Recon"`

`<params>` is

    { "Tree":<tree>,
      "Alignment":<alignment>,
      <optional-args> }

or 

    { "POGTree":<POGtree>, 
      <optional-args>}

and

`<optional-args>` is

    "Inference":"Joint", (or) "Inference":"Marginal","Ancestor":<anc-ID>,
    "Indels":<indel-method>,
    "Model":<subst-model>,
    "Rates":[<rate-pos1>,<rate-pos2>, ...]

where `<indel-method>` is `"BEP"` (default), `"BEML"`, `"SICP"`, `"SICML"`, `"PSP"`, or `"PSML"`;
`<subst-model>` is `"JTT"` (default), `"Dayhoff"`, `"LG"`, `"WAG"`, `"JC"`, or `"Yang"`.

Result is
    
    "GRASP_version":<GRASP-version>,
    "Ancestors":[<list-of-POGs>]

#### Example

    { "Command":"Recon",
      "Auth":"Guest",
      "Params":
        { "Alignment":
            { "Sequences":[{"Seq":[null,null,null,"L",null],"Name":"sequence38"},{"Seq":[null,null,null,null,"R"],"Name":"sequence64"},{"Seq":[null,null,"S",null,"R"],"Name":"sequence87"},{"Seq":[null,null,"A",null,"C"],"Name":"sequence77"},{"Seq":[null,null,"A",null,"C"],"Name":"sequence110"},{"Seq":[null,null,null,"V","T"],"Name":"sequence239"},{"Seq":[null,null,null,"I",null],"Name":"sequence101"},{"Seq":[null,null,"A",null,"E"],"Name":"sequence203"},{"Seq":[null,null,"A",null,"H"],"Name":"sequence30"},{"Seq":[null,null,"R",null,"H"],"Name":"sequence99"},{"Seq":[null,null,"A",null,"H"],"Name":"sequence60"},{"Seq":[null,"T","T",null,"F"],"Name":"sequence50"},{"Seq":[null,"T","T",null,"F"],"Name":"sequence111"},{"Seq":[null,"T","A",null,"Y"],"Name":"sequence197"},{"Seq":[null,null,"A",null,"I"],"Name":"sequence72"},{"Seq":[null,null,"A",null,"N"],"Name":"sequence229"},{"Seq":[null,null,null,"V","D"],"Name":"sequence3"},{"Seq":[null,null,null,"I","D"],"Name":"sequence125"},{"Seq":[null,null,"L",null,"H"],"Name":"sequence86"},{"Seq":[null,null,"S",null,"H"],"Name":"sequence237"},{"Seq":[null,null,"A",null,"H"],"Name":"sequence79"},{"Seq":[null,null,"A",null,"H"],"Name":"sequence75"},{"Seq":[null,null,"A",null,"H"],"Name":"sequence62"},{"Seq":[null,null,"A",null,"H"],"Name":"sequence209"},{"Seq":[null,null,"A",null,"R"],"Name":"sequence117"},{"Seq":[null,null,"A",null,"K"],"Name":"sequence191"},{"Seq":[null,null,"A",null,"Q"],"Name":"sequence56"},{"Seq":[null,null,"A",null,"H"],"Name":"sequence137"},{"Seq":[null,null,"A",null,"H"],"Name":"sequence1"},{"Seq":[null,null,"A",null,"H"],"Name":"sequence66"},{"Seq":[null,null,"Q",null,"Q"],"Name":"sequence212"},{"Seq":["L","G","E",null,"H"],"Name":"sequence108"},{"Seq":["L","G","E",null,"H"],"Name":"sequence228"},{"Seq":[null,"T","E",null,"F"],"Name":"sequence68"},{"Seq":[null,"T","A",null,"F"],"Name":"sequence216"},{"Seq":[null,null,"A",null,"C"],"Name":"sequence95"},{"Seq":[null,null,"A",null,"R"],"Name":"sequence113"},{"Seq":[null,null,"A",null,"R"],"Name":"sequence73"},{"Seq":[null,null,"A",null,"R"],"Name":"sequence57"},{"Seq":[null,null,"A",null,"D"],"Name":"sequence6"},{"Seq":[null,null,"A",null,"D"],"Name":"sequence52"},{"Seq":[null,null,"T",null,"H"],"Name":"sequence139"},{"Seq":[null,null,"E",null,"H"],"Name":"sequence47"},{"Seq":[null,null,"E",null,"H"],"Name":"sequence234"},{"Seq":[null,null,"A",null,"T"],"Name":"sequence159"},{"Seq":[null,null,"A",null,"C"],"Name":"sequence247"},{"Seq":[null,null,"A",null,"N"],"Name":"sequence39"},{"Seq":[null,null,"A",null,"N"],"Name":"sequence8"},{"Seq":[null,null,"A",null,"N"],"Name":"sequence119"},{"Seq":[null,null,"A",null,"R"],"Name":"sequence33"},{"Seq":[null,null,"A",null,"R"],"Name":"sequence210"},{"Seq":[null,null,"A",null,"H"],"Name":"sequence43"},{"Seq":[null,null,"E",null,"K"],"Name":"sequence121"},{"Seq":[null,null,"E",null,"K"],"Name":"sequence92"},{"Seq":[null,null,"S",null,"K"],"Name":"sequence27"},{"Seq":[null,null,"A",null,"K"],"Name":"sequence170"},{"Seq":[null,null,"G",null,"H"],"Name":"sequence231"},{"Seq":[null,null,"A",null,"H"],"Name":"sequence54"},{"Seq":[null,null,"A",null,"N"],"Name":"sequence168"},{"Seq":[null,null,"A",null,"N"],"Name":"sequence219"},{"Seq":[null,null,"T",null,"R"],"Name":"sequence7"},{"Seq":[null,null,"T",null,"R"],"Name":"sequence199"},{"Seq":[null,null,"T",null,"R"],"Name":"sequence97"},{"Seq":[null,null,"T",null,"R"],"Name":"sequence74"},{"Seq":[null,null,"A",null,"H"],"Name":"sequence134"},{"Seq":[null,null,"A",null,"H"],"Name":"sequence179"}],
              "Datatype":{"Predef":"Protein"} },
          "Tree":
            { "Parents":[-1,0,1,2,3,4,5,6,7,7,6,5,11,11,4,3,15,16,16,15,19,19,21,21,2,24,25,26,27,28,28,27,26,25,33,33,35,36,37,37,36,35,41,41,24,44,44,1,47,48,49,50,50,49,48,54,54,47,0,58,59,59,61,62,63,64,64,66,66,63,69,69,62,72,72,74,75,76,76,75,79,79,74,82,83,83,82,61,87,88,88,90,90,87,93,93,58,96,97,98,98,100,100,97,103,104,105,106,106,105,109,109,104,103,113,113,96,116,117,117,119,119,116,122,123,124,125,125,124,123,122],
              "Labels":["0","1","2","3","4","5","6","7","sequence47","sequence234","sequence212","8","sequence108","sequence228","sequence87","9","10","sequence50","sequence111","11","sequence197","12","sequence68","sequence216","13","14","15","16","17","sequence159","sequence247","sequence95","sequence72","18","sequence229","19","20","21","sequence168","sequence219","sequence39","22","sequence8","sequence119","23","sequence77","sequence110","24","25","26","27","sequence3","sequence125","sequence239","28","sequence101","sequence203","sequence38","29","30","sequence64","31","32","33","34","sequence113","35","sequence33","sequence210","36","sequence73","sequence57","37","sequence86","38","39","40","sequence7","sequence199","41","sequence97","sequence74","42","43","sequence134","sequence179","sequence43","44","45","sequence237","46","sequence6","sequence52","47","sequence79","sequence75","48","49","50","sequence30","51","sequence62","sequence209","52","53","54","55","sequence121","sequence92","56","sequence27","sequence170","sequence117","57","sequence191","sequence56","58","59","sequence99","60","sequence137","sequence1","61","62","63","64","sequence231","sequence54","sequence139","sequence66","sequence60"],
              "Distances":[0,0.086838,0.023163,0.026222,0.342957,0.02521,0.051304,0.174255,0.145796,0.145796,0.320051,0.284419,0.086936,0.086936,0.396566,0.477527,0.258728,0.003268,0.003268,0.142402,0.119593,0.040549,0.079044,0.079044,0.169549,0.077884,0.066855,0.092639,0.041577,0.317241,0.317241,0.358819,0.451456,0.392445,0.125867,0.016502,0.05269,0.054446,0.002228,0.002228,0.056674,0.101538,0.007827,0.007827,0.361893,0.234304,0.234304,0.050065,0.02764,0.248498,0.3683,0.094405,0.094405,0.462704,0.048629,0.662573,0.662573,0.738842,0.058064,0.0922,0.725482,0.077351,0.028822,0.270835,0.164158,0.184316,0.121621,0.062695,0.062695,0.234931,0.113543,0.113543,0.08903,0.530278,0.079736,0.262731,0.11578,0.072031,0.072031,0.031693,0.156119,0.156119,0.352005,0.064503,0.034035,0.034035,0.098538,0.192432,0.014164,0.441535,0.44106,4.74E-4,4.74E-4,0.246016,0.209683,0.209683,0.081914,0.014711,0.473628,0.247428,0.019549,0.22788,0.22788,0.152647,0.180719,0.005074,0.365646,0.016972,0.016972,0.122938,0.25968,0.25968,0.387692,0.119619,0.448791,0.448791,0.240317,0.012407,0.483044,0.099375,0.38367,0.38367,0.071027,0.021988,0.161969,0.129501,0.110966,0.110966,0.240467,0.402436,0.424423],
              "Branchpoints":131 } } }

which results in 

    { "Job":14,
      "Result":
        { "GRASP_version":"12-Dec-2022",
          "Ancestors":[
            { "Adjacent":[[4],[]],"Starts":[2],"Edgeindices":[[2,4],[-1,2],[4,5]],"Nodetype":"class dat.pog.SymNode","Size":5,"Indices":[2,4],"GRASP_version":"12-Dec-2022","Directed":true,"Terminated":true,"Edges":[{"Recip":true,"Backward":true,"Forward":true,"Weight":0},{"Recip":true,"Backward":true,"Forward":true,"Weight":0},{"Recip":true,"Backward":true,"Forward":true,"Weight":0}],"Nodes":[{"Value":"A"},{"Value":"H"}],"Datatype":"class dat.pog.POGraph","Name":"0","Ends":[4],"Edgetype":"class dat.pog.POGraph$BidirEdge" },
            { "Adjacent":[[4],[]],"Starts":[2],"Edgeindices":[[2,4],[-1,2],[4,5]],"Nodetype":"class dat.pog.SymNode","Size":5,"Indices":[2,4],"GRASP_version":"12-Dec-2022","Directed":true,"Terminated":true,"Edges":[{"Recip":true,"Backward":true,"Forward":true,"Weight":0},{"Recip":true,"Backward":true,"Forward":true,"Weight":0},{"Recip":true,"Backward":true,"Forward":true,"Weight":0}],"Nodes":[{"Value":"A"},{"Value":"H"}],"Datatype":"class dat.pog.POGraph","Name":"1","Ends":[4],"Edgetype":"class dat.pog.POGraph$BidirEdge" },
            { "Adjacent":[[4],[]],"Starts":[2],"Edgeindices":[[2,4],[-1,2],[4,5]],"Nodetype":"class dat.pog.SymNode","Size":5,"Indices":[2,4],"GRASP_version":"12-Dec-2022","Directed":true,"Terminated":true,"Edges":[{"Recip":true,"Backward":true,"Forward":true,"Weight":0},{"Recip":true,"Backward":true,"Forward":true,"Weight":0},{"Recip":true,"Backward":true,"Forward":true,"Weight":0}],"Nodes":[{"Value":"A"},{"Value":"H"}],"Datatype":"class dat.pog.POGraph","Name":"2","Ends":[4],"Edgetype":"class dat.pog.POGraph$BidirEdge" },
            { ... } ... ] } }

Note: consider returning the full package with input data, i.e.

    { "Prediction":{"Input":{<input-JSON>},
                    "Ancestors":[<list-of-POGs>],
                    "Datatype":"Prediction"}}}

### Train BP2Prop model: `<command>` is `"Train"`

`<params>` is

    { "Tree":<tree>,
      "Dataset":<dataset>,
      "States":[<state1>,<state2>,...],
      <optional-args> }

where `<dataset>` specifies a data set (see below).

and

`<optional-args>` is

      "Distrib":<distrib>,
      "Leaves-only":<true/false>,
      "Rate":<rate>,
      "Seed":<seed>,
      "Gamma":<gamma>

where `<distrib>` specifies a distribution appropriate for the values in the `<dataset>` (see below).

Result is a new/refined distribution (see example below).

#### Example

    { "Command":"Train",
      "Auth":"Guest",
      "Params":
        { "States":["A","B"],
          "Dataset":
            { "Headers":["S009","S005","S002","S006","S003","S001","S008","S010","S004","S007"],
              "Data":[[3.63],[3.81],[2.89],[3.81],[2.54],[2.76],[3.79],[3.7],[1.94],[3.97]]},
              "Tree":<tree> } } }

which results in
    
    { "Distrib":
        { "Condition":[["A"],["B"]],
          "Pr":[[3.784926135903969,0.056738891699391655],[2.5324588293595744,0.056738891699391655]],
          "Index":[0,1],
          "Domain":"dat.Continuous@3bd5adde" } }

### Infer properties: `<command>` is `"Infer"`

`<params>` is

    { "Tree":<tree>,
      "Dataset":<dataset>,
      "States":[<state1>,<state2>,...],
      "Distrib":<distrib>,
      <optional-args> }

and

`<optional-args>` is

      "Distrib":<distrib>,
      "Leaves-only":<true/false>,
      "Rate":<rate>,
      "Seed":<seed>,
      "Gamma":<gamma>

where `<distrib>` specifies a distribution appropriate for the values in the `<dataset>` (see below).

Result is a new/refined distribution (see example below).

#### Example

    { "Command":"Infer",
      "Auth":"Guest",
      "Params":
        { "States":["A","B"],
          "Leaves-only":true,
          "Dataset":
            { "Headers":["S009","S005","S002","S006","S003","S001","S008","S010","S004","S007"],"Data":[[3.63,3.33],[3.81,3.21],[2.89,2.93],[3.81,3.51],[2.54,2.59],[2.76,2.96],[3.79,3.49],[3.7,3.4],[1.94,2.24],[3.97,3.44]] },
          "Tree":<tree>,
          "Inference":"Marginal",
          "Distrib":
            { "Condition":[["A"],["B"]],"Pr":[[3.784926135903969,0.056738891699391655],[2.532458829359575,0.056738891699391655]],"Index":[0,1],"Domain":"dat.Continuous@3cf71bc7"},
          "Ancestor":0 } }

Server responds:

    { "N0":[
        { "Pr":[0.6652270145537978,0.3347729854462022],"Domain":{"Size":2,"Values":["A","B"],"Datatype":"String"} },
        { "Pr":[0.649968113685095,0.350031886314905],"Domain":{"Size":2,"Values":["A","B"],"Datatype":"String"} }] }

## Data structures

### `<tree>` (`dat.phylo.IdxTree`)

    { "Branchpoints":<number-of-idxs>,
      "Labels":[<list-of-labels>],
      "Parents":[<list-of-parent-idxs>],     // -1 represents no parent
      "Distances":[<list-of-dists>] }

#### Example

Newick definition

    ((A:0.6,((B:3.3,(C:1.0,D:2.5)cd:1.8)bcd:5,((E:3.9,F:4.5)ef:2.5,G:0.3)efg:7)X:3.2)Y:0.5,H:1.1)I:0.2

JSON definition (note internal nodes are re-numbered, left-to-right, depth-first)

    {"Parents":[-1,0,1,1,3,4,4,6,6,3,9,10,10,9,0],"Labels":["0","1","A","2","3","B","4","C","D","5","6","E","F","G","H"],"Distances":[0,0.5,0.6,3.2,5,3.3,1.8,1,2.5,7,2.5,3.9,4.5,0.3,1.1],"Branchpoints":15}

### `<alignment>` (`dat.EnumSeq.Alignment`)


### `<POG>` (`dat.pog.POGraph`)

#### Example

Graph with a maximum of 5 nodes, here with two (indices 2 and 4), where 2 is adjacent to 4. 
The graph starts with node with index 2 and ends with node with index 4.

    { "Size":5,
      "Indices":[2,4],
      "Adjacent":[[4],[]],
      "Starts":[2],
      "Ends":[4],
      "Nodetype":"class dat.pog.SymNode",
      "GRASP_version":"12-Dec-2022",
      "Directed":true,
      "Terminated":true,
      "Nodes":[{"Value":"E"},{"Value":"H"}],
      "Datatype":"class dat.pog.POGraph",
      "Name":"sequence47" }

Additional fields in graphs with edge attributes (including POGs):

      "Edgeindices":[[2,4],[-1,2],[4,5]],
      "Edges":[{"Recip":true,"Backward":true,"Forward":true,"Weight":0}, ...]

### `<POGTree>` (`dat.pog.POGTree`)

    { "Tree":<tree>,
      "Hashcode":-631886170,
      "Extants":[<POG1>, <POG2>, ...] }

### `<dataset>` (`api.JSONUtils.DataSet`)

    { "Headers":[<header1>,<header2>, ...],
      "Data":[ [ ... ] ] }

Note: the data matrix is given column 1, column 2, etc, with values in row 1, row 2 etc within

### `<bnode>` (`bn.BNode` i.e. Bayesian network node) 

Here's a discrete/enumerable CPT node (with one enumerable parent variable):

    { "Condition":[["A"],["B"],["C"]],
      "Pr":[[0.36,0.38,0.26],[0.50,0.32,0.18],[0.07,0.39,0.54]],
      "Index":[0,1,2],
      "Domain":
        { "Size":3,
          "Values":[0,1,2],
          "Datatype":"Integer" } }

Here's a continuous/Gaussian GDT node (with the same enumerable parent variable):

    { "Condition":[["A"],["B"],["C"]],
      "Pr":[[2.1,0.5],[6.4,0.5],[1.3,0.5]],
      "Index":[0,1,2],
      "Domain":"dat.Continuous@525b461a" }











