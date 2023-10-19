# JSON API 

## Requests for `asr.GServer`

    { "Command":<request>,
      "Job":<job-number> }

### `<request>` is `"Retrieve"`
Recover the original request under the specified job-number

### `<request>` is `"Output"`
Return the result of the request under the specified job-number
on the format

    { "Job":<job-number>,
      "Result":{<result-JSON>} }

### `<request>` is `"Place"`
Find the place in queue of the request under the specified job-number

### List all jobs

    { "Command":"Status"}

returns

    { "Jobs":[<list-of-jobs>],
      "Clients":<number-of-clients>}

for example

    { "Jobs":[
        { "Status":"COMPLETED","Threads":1,"Command":"Fake","Priority":0,"Memory":1,"Auth":"Guest","Job":1,"Place":0 },
        { "Status":"COMPLETED","Threads":1,"Command":"Fake","Priority":0,"Memory":1,"Auth":"Guest","Job":2,"Place":0 },
        { "Status":"RUNNING","Threads":1,"Command":"Fake","Priority":0,"Memory":1,"Auth":"Guest","Job":3,"Place":0 },
        { "Status":"WAITING","Threads":1,"Command":"Fake","Priority":0,"Memory":1,"Auth":"Guest","Job":4,"Place":1 }],
      "Clients":1 }

## Commands for `asr.GServer`
Execute/submit job; a "command" is a request that will require the server to dedicate resources both in terms of memory and CPU time, 
which is why it is allocated a "job-number" and typically queued 

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

As a one-line JSON-string:

    {"Command":"Recon","Auth":"Guest","Params":{"Alignment":{"Sequences":[{"Seq":[null,null,null,"L",null],"Name":"sequence38"},{"Seq":[null,null,null,null,"R"],"Name":"sequence64"},{"Seq":[null,null,"S",null,"R"],"Name":"sequence87"},{"Seq":[null,null,"A",null,"C"],"Name":"sequence77"},{"Seq":[null,null,"A",null,"C"],"Name":"sequence110"},{"Seq":[null,null,null,"V","T"],"Name":"sequence239"},{"Seq":[null,null,null,"I",null],"Name":"sequence101"},{"Seq":[null,null,"A",null,"E"],"Name":"sequence203"},{"Seq":[null,null,"A",null,"H"],"Name":"sequence30"},{"Seq":[null,null,"R",null,"H"],"Name":"sequence99"},{"Seq":[null,null,"A",null,"H"],"Name":"sequence60"},{"Seq":[null,"T","T",null,"F"],"Name":"sequence50"},{"Seq":[null,"T","T",null,"F"],"Name":"sequence111"},{"Seq":[null,"T","A",null,"Y"],"Name":"sequence197"},{"Seq":[null,null,"A",null,"I"],"Name":"sequence72"},{"Seq":[null,null,"A",null,"N"],"Name":"sequence229"},{"Seq":[null,null,null,"V","D"],"Name":"sequence3"},{"Seq":[null,null,null,"I","D"],"Name":"sequence125"},{"Seq":[null,null,"L",null,"H"],"Name":"sequence86"},{"Seq":[null,null,"S",null,"H"],"Name":"sequence237"},{"Seq":[null,null,"A",null,"H"],"Name":"sequence79"},{"Seq":[null,null,"A",null,"H"],"Name":"sequence75"},{"Seq":[null,null,"A",null,"H"],"Name":"sequence62"},{"Seq":[null,null,"A",null,"H"],"Name":"sequence209"},{"Seq":[null,null,"A",null,"R"],"Name":"sequence117"},{"Seq":[null,null,"A",null,"K"],"Name":"sequence191"},{"Seq":[null,null,"A",null,"Q"],"Name":"sequence56"},{"Seq":[null,null,"A",null,"H"],"Name":"sequence137"},{"Seq":[null,null,"A",null,"H"],"Name":"sequence1"},{"Seq":[null,null,"A",null,"H"],"Name":"sequence66"},{"Seq":[null,null,"Q",null,"Q"],"Name":"sequence212"},{"Seq":["L","G","E",null,"H"],"Name":"sequence108"},{"Seq":["L","G","E",null,"H"],"Name":"sequence228"},{"Seq":[null,"T","E",null,"F"],"Name":"sequence68"},{"Seq":[null,"T","A",null,"F"],"Name":"sequence216"},{"Seq":[null,null,"A",null,"C"],"Name":"sequence95"},{"Seq":[null,null,"A",null,"R"],"Name":"sequence113"},{"Seq":[null,null,"A",null,"R"],"Name":"sequence73"},{"Seq":[null,null,"A",null,"R"],"Name":"sequence57"},{"Seq":[null,null,"A",null,"D"],"Name":"sequence6"},{"Seq":[null,null,"A",null,"D"],"Name":"sequence52"},{"Seq":[null,null,"T",null,"H"],"Name":"sequence139"},{"Seq":[null,null,"E",null,"H"],"Name":"sequence47"},{"Seq":[null,null,"E",null,"H"],"Name":"sequence234"},{"Seq":[null,null,"A",null,"T"],"Name":"sequence159"},{"Seq":[null,null,"A",null,"C"],"Name":"sequence247"},{"Seq":[null,null,"A",null,"N"],"Name":"sequence39"},{"Seq":[null,null,"A",null,"N"],"Name":"sequence8"},{"Seq":[null,null,"A",null,"N"],"Name":"sequence119"},{"Seq":[null,null,"A",null,"R"],"Name":"sequence33"},{"Seq":[null,null,"A",null,"R"],"Name":"sequence210"},{"Seq":[null,null,"A",null,"H"],"Name":"sequence43"},{"Seq":[null,null,"E",null,"K"],"Name":"sequence121"},{"Seq":[null,null,"E",null,"K"],"Name":"sequence92"},{"Seq":[null,null,"S",null,"K"],"Name":"sequence27"},{"Seq":[null,null,"A",null,"K"],"Name":"sequence170"},{"Seq":[null,null,"G",null,"H"],"Name":"sequence231"},{"Seq":[null,null,"A",null,"H"],"Name":"sequence54"},{"Seq":[null,null,"A",null,"N"],"Name":"sequence168"},{"Seq":[null,null,"A",null,"N"],"Name":"sequence219"},{"Seq":[null,null,"T",null,"R"],"Name":"sequence7"},{"Seq":[null,null,"T",null,"R"],"Name":"sequence199"},{"Seq":[null,null,"T",null,"R"],"Name":"sequence97"},{"Seq":[null,null,"T",null,"R"],"Name":"sequence74"},{"Seq":[null,null,"A",null,"H"],"Name":"sequence134"},{"Seq":[null,null,"A",null,"H"],"Name":"sequence179"}],"Datatype":{"Predef":"Protein"}},"Tree":{"Parents":[-1,0,1,2,3,4,5,6,7,7,6,5,11,11,4,3,15,16,16,15,19,19,21,21,2,24,25,26,27,28,28,27,26,25,33,33,35,36,37,37,36,35,41,41,24,44,44,1,47,48,49,50,50,49,48,54,54,47,0,58,59,59,61,62,63,64,64,66,66,63,69,69,62,72,72,74,75,76,76,75,79,79,74,82,83,83,82,61,87,88,88,90,90,87,93,93,58,96,97,98,98,100,100,97,103,104,105,106,106,105,109,109,104,103,113,113,96,116,117,117,119,119,116,122,123,124,125,125,124,123,122],"Labels":["0","1","2","3","4","5","6","7","sequence47","sequence234","sequence212","8","sequence108","sequence228","sequence87","9","10","sequence50","sequence111","11","sequence197","12","sequence68","sequence216","13","14","15","16","17","sequence159","sequence247","sequence95","sequence72","18","sequence229","19","20","21","sequence168","sequence219","sequence39","22","sequence8","sequence119","23","sequence77","sequence110","24","25","26","27","sequence3","sequence125","sequence239","28","sequence101","sequence203","sequence38","29","30","sequence64","31","32","33","34","sequence113","35","sequence33","sequence210","36","sequence73","sequence57","37","sequence86","38","39","40","sequence7","sequence199","41","sequence97","sequence74","42","43","sequence134","sequence179","sequence43","44","45","sequence237","46","sequence6","sequence52","47","sequence79","sequence75","48","49","50","sequence30","51","sequence62","sequence209","52","53","54","55","sequence121","sequence92","56","sequence27","sequence170","sequence117","57","sequence191","sequence56","58","59","sequence99","60","sequence137","sequence1","61","62","63","64","sequence231","sequence54","sequence139","sequence66","sequence60"],"Distances":[0,0.086838,0.023163,0.026222,0.342957,0.02521,0.051304,0.174255,0.145796,0.145796,0.320051,0.284419,0.086936,0.086936,0.396566,0.477527,0.258728,0.003268,0.003268,0.142402,0.119593,0.040549,0.079044,0.079044,0.169549,0.077884,0.066855,0.092639,0.041577,0.317241,0.317241,0.358819,0.451456,0.392445,0.125867,0.016502,0.05269,0.054446,0.002228,0.002228,0.056674,0.101538,0.007827,0.007827,0.361893,0.234304,0.234304,0.050065,0.02764,0.248498,0.3683,0.094405,0.094405,0.462704,0.048629,0.662573,0.662573,0.738842,0.058064,0.0922,0.725482,0.077351,0.028822,0.270835,0.164158,0.184316,0.121621,0.062695,0.062695,0.234931,0.113543,0.113543,0.08903,0.530278,0.079736,0.262731,0.11578,0.072031,0.072031,0.031693,0.156119,0.156119,0.352005,0.064503,0.034035,0.034035,0.098538,0.192432,0.014164,0.441535,0.44106,4.74E-4,4.74E-4,0.246016,0.209683,0.209683,0.081914,0.014711,0.473628,0.247428,0.019549,0.22788,0.22788,0.152647,0.180719,0.005074,0.365646,0.016972,0.016972,0.122938,0.25968,0.25968,0.387692,0.119619,0.448791,0.448791,0.240317,0.012407,0.483044,0.099375,0.38367,0.38367,0.071027,0.021988,0.161969,0.129501,0.110966,0.110966,0.240467,0.402436,0.424423],"Branchpoints":131}}}

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

      "Distrib":<bnode>,
      "Leaves-only":<true/false>,
      "Rate":<rate>,
      "Seed":<seed>,
      "Gamma":<gamma>

where `<bnode>` specifies a probability distribution appropriate for the values in the `<dataset>` 
(see below, `bnode` is actually a part of a Bayesian network). 
If you specify one, it will be used as a starting point for training.

The result is a new/refined distribution (see example below).

#### Example (continuous)

    { "Command":"Train",
      "Auth":"Guest",
      "Params":
        { "States":["A","B"],
          "Dataset":
            { "Headers":["S009","S005","S002","S006","S003","S001","S008","S010","S004","S007"],
              "Data":[[3.63],[3.81],[2.89],[3.81],[2.54],[2.76],[3.79],[3.7],[1.94],[3.97]]},
              "Tree":<tree> } } }

which results in a Gaussian distribution, or `bnode` (mean and variance for each component, or latent value): 
    
    { "Distrib":
        { "Condition":[["A"],["B"]],
          "Pr":[[3.784926135903969,0.056738891699391655],[2.5324588293595744,0.056738891699391655]],
          "Index":[0,1],
          "Domain":"dat.Continuous@3bd5adde" } }

#### Example (discrete)

    { "Command":"Train",
      "Auth":"Guest",
      "Params":
        { "States":["A","B","C"],
          "Dataset":
            { "Headers":["H9B4I9","B8Y1H0","A8VWC5","Q47R94","Q1KLC8","P29957","B6RB08","P04746","P00690","H2N0D4","A0SEG1","Q2KJQ1","B8Y698","P08144","D4P4Y7","T1WDH2","T1WE96","Q4A3E0","Q8LJQ6","P00693","P17654","P04063","O33476","O08452","Q2QC88","O93647","Q6WUB6","A0A060DAC6","B1VK33","P06279","P00692","P06278","Q5UZY3","L8B068","D8J7H2","A5ILB0","P96107","Q8A1G3","P20845"],
              "Data":[["G"],["G"],["G"],["Q"],["Q"],["Q"],["Q"],["Q"],["Q"],["Q"],["Q"],["Q"],["H"],["Q"],["G"],["G"],["G"],["G"],["A"],["A"],["A"],["A"],["G"],["G"],["G"],["G"],["A"],["G"],["G"],["G"],["A"],["A"],["L"],["D"],["G"],["G"],["G"],["G"],["G"]] },
          "Tree":
            { "Parents":[-1,0,0,0,3,3,5,6,6,8,9,9,8,12,12,14,15,16,17,17,16,20,20,22,22,24,24,26,26,15,29,30,30,32,32,29,35,35,37,37,14,40,41,41,43,43,45,45,47,47,49,49,40,52,53,54,54,53,57,57,52,60,60,62,62,64,64,5,67,67,69,70,70,69,73,73],
              "Labels":["0","P08144","B8Y698","1","Q2KJQ1","2","3","P29957","4","5","Q47R94","Q1KLC8","6","Q4A3E0","7","8","9","10","Q5UZY3","L8B068","11","Q8A1G3","12","D8J7H2","13","P20845","14","A5ILB0","P96107","15","16","D4P4Y7","17","T1WDH2","T1WE96","18","H9B4I9","19","B8Y1H0","A8VWC5","20","21","Q6WUB6","22","A0A060DAC6","23","B1VK33","24","P06279","25","P00692","P06278","26","27","28","O33476","O08452","29","Q2QC88","O93647","30","Q8LJQ6","31","P00693","32","P17654","P04063","33","B6RB08","34","35","P04746","P00690","36","H2N0D4","A0SEG1"],
              "Distances":[0,0.243543659,0.328253867,0.061234395,0.302246327,0.104112732,0.129975033,0.492317028,0.168477423,0.472234961,5.0E-9,0.001470602,0.838140064,0.758927261,0.118993954,0.217610635,0.256851036,0.71053731,0.388091681,0.525286732,0.279916536,0.876914346,0.093326396,0.897751479,0.091555295,0.581536564,0.645218372,0.009417249,0.013626135,0.388197737,0.242001003,0.554412114,0.442555732,0.198231144,0.203377051,1.173898406,0.051873954,0.030364499,0.002987842,0.002295265,0.692720871,0.165853581,0.630412028,0.301179158,0.564737357,0.098304443,0.525139427,0.114171367,0.243352437,0.187788353,0.148531194,0.088716021,0.144788807,0.522130083,0.035341913,0.052453145,0.084658893,0.032451459,0.089455947,0.071707337,1.039127746,0.133051046,0.105731012,0.187271089,0.078869174,0.093152494,0.082644055,0.052339221,0.342080258,0.20996685,0.115876665,0.090541678,0.04644842,0.068089852,0.0917515,0.142809732],
              "Branchpoints":76 } } }

which results in a discrete distribution (multinomial distributions, one for each component)

    { "Distrib":
        { "Condition":[["A"],["B"],["C"]],
          "Pr":[
            [0.91,      5.98E-79, 2.13E-27, 1.08E-39, 0.09,     3.29E-28],
            [2.80E-122, 8.16E-65, 4.59E-31, 1.00,     1.76E-97, 2.65E-28],
            [3.53E-91,  0.78,     0.11,     2.13E-13, 3.65E-65, 0.11    ] ],
          "Index":[0,1,2],
          "Domain":
            { "Size":6,
              "Values":["Q","A","D","G","H","L"],
              "Datatype":"String" } } }

### Infer properties: `<command>` is `"Infer"`

`<params>` is

    { "Tree":<tree>,
      "Dataset":<dataset>,
      "States":[<state1>,<state2>,...],
      "Distrib":<bnode>,
      <optional-args> }

where `<bnode>` must specify a distribution appropriate for the values in the `<dataset>` (see below); 
`<optional-args>` is

      "Leaves-only":<true/false>,
      "Rate":<rate>,
      "Seed":<seed>,
      "Gamma":<gamma>

Result is a new/refined distribution (see example below).

#### Example of marginal reconstruction of discrete states

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

#### Example with joint reconstruction of discrete states

    { "Command":"Infer",
      "Auth":"Guest",
      "Params":
        { "States":["a","b","c","d"],
          "Leaves-only":true,
          "Dataset":
            { "Headers":["H9B4I9","B8Y1H0","A8VWC5","Q47R94","Q1KLC8","P29957","B6RB08","P04746","P00690","H2N0D4","A0SEG1","Q2KJQ1","B8Y698","P08144","D4P4Y7","T1WDH2","T1WE96","Q4A3E0","Q8LJQ6","P00693","P17654","P04063","O33476","O08452","Q2QC88","O93647","Q6WUB6","A0A060DAC6","B1VK33","P06279","P00692","P06278","Q5UZY3","L8B068","D8J7H2","A5ILB0","P96107","Q8A1G3","P20845"],
              "Data":[["G"],["G"],["G"],["Q"],["Q"],["Q"],["Q"],["Q"],["Q"],["Q"],["Q"],["Q"],["H"],["Q"],["G"],["G"],["G"],["G"],["A"],["A"],["A"],["A"],["G"],["G"],["G"],["G"],["A"],["G"],["G"],["G"],["A"],["A"],["L"],["D"],["G"],["G"],["G"],["G"],["G"] ] },
          "Tree":{<tree>},
          "Inference":"Joint",
          "Distrib":{"Condition":[["a"],["b"],["c"],["d"]],"Pr":[[2.28E-71,1.71E-44,0.12,0.75,5.72E-55,0.12],[8.51E-84,5.14E-41,2.35E-34,1,7.81E-67,2.14E-32],[2.30E-70,0.99,4.48E-28,1.40E-11,1.78E-52,2.87E-28],[0.91,8.09E-66,4.53E-36,1.96E-42,0.09,4.35E-37]],"Index":[0,1,2,3],"Domain":{"Size":6,"Values":["Q","A","D","G","H","L"],"Datatype":"String"} } } }

Server responds with a (completed) dataset:

    { "Predict":
        { "Headers":["N0","P08144","B8Y698","N1","Q2KJQ1","N2","N3","P29957","N4","N5","Q47R94","Q1KLC8","N6","Q4A3E0","N7","N8","N9","N10","Q5UZY3","L8B068","N11","Q8A1G3","N12","D8J7H2","N13","P20845","N14","A5ILB0","P96107","N15","N16","D4P4Y7","N17","T1WDH2","T1WE96","N18","H9B4I9","N19","B8Y1H0","A8VWC5","N20","N21","Q6WUB6","N22","A0A060DAC6","N23","B1VK33","N24","P06279","N25","P00692","P06278","N26","N27","N28","O33476","O08452","N29","Q2QC88","O93647","N30","Q8LJQ6","N31","P00693","N32","P17654","P04063","N33","B6RB08","N34","N35","P04746","P00690","N36","H2N0D4","A0SEG1"],
          "Data":[["d"],["Q"],["H"],["d"],["Q"],["d"],["d"],["Q"],["d"],["d"],["Q"],["Q"],["b"],["G"],["b"],["b"],["b"],["a"],["L"],["D"],["b"],["G"],["b"],["G"],["b"],["G"],["b"],["G"],["G"],["b"],["b"],["G"],["b"],["G"],["G"],["b"],["G"],["b"],["G"],["G"],["b"],["b"],["A"],["b"],["G"],["b"],["G"],["b"],["G"],["c"],["A"],["A"],["b"],["b"],["b"],["G"],["G"],["b"],["G"],["G"],["c"],["A"],["c"],["A"],["c"],["A"],["A"],["d"],["Q"],["d"],["d"],["Q"],["Q"],["d"],["Q"],["Q"]] } }

### Train "modes" (example)

    { "Command":"TrainModes",
      "Auth":"Guest",
      "Params":
        { "Gamma":1,
          "Rounds":10,
          "Dataset":
            { "Items":["P25910","Q704V1",...],
              "Features":["Pos184","Pos186","Pos306"],
              "Data":[[["H","H","H"],["H","H","H"], ...] ] },
          "Seed":3,
          "Rate":1,
          "Tree":
            { "Parents":[-1,0, ...],
              "Labels":["0","P25910", ...],
              "Distances":[0,0.73, ...],
              "Branchpoints":220 },
          "Distrib":
            { "Targets":[[0],[0],[0]],
              "Modetypes":[{"Size":3,"Values":["A","B","C"],"Datatype":"Character"}],
              "Nodes":[
                { "Condition":[],
                  "Pr":[],
                  "Variable":{"Domain":{"Predef":"Protein"},"Name":"N0__Pos184"},
                  "Nodetype":"CPT",
                  "Index":[] },
                { "Condition":[],
                  "Pr":[],
                  "Variable":{"Domain":{"Predef":"Protein"},"Name":"N0__Pos186"},
                  "Nodetype":"CPT",
                  "Index":[] },
                { "Condition":[],
                  "Pr":[],
                  "Variable":{"Domain":{"Predef":"Protein"},"Name":"N0__Pos306"},
                  "Nodetype":"CPT",
                  "Index":[] } ],
              "Name":"N0" } } }

### Infer "modes" (example)

    { "Command":"InferModes",
      "Auth":"Guest",
      "Params":
        { "Gamma":1,
          "Latent":true,
          "Rounds":10,
          "Leaves-only":false,
          "Dataset":
            { "Items":["P25910","Q704V1", ...],
              "Features":["Pos184","Pos186","Pos306"],
              "Data":[[["H","H","H"],["H","H","H"], ...] ] },
          "Seed":3,
          "Rate":1,
          "Tree":
            { "Parents":[-1,0, ...],
              "Labels":["0","P25910", ...],
              "Distances":[0,0.73, ...],
              "Branchpoints":220 },
          "Inference":"Marginal",
          "Distrib":
            { "Targets":[[0],[0],[0]],
              "Modetypes":[{"Size":3,"Values":["A","B","C"],"Datatype":"String"}],
              "Nodes":[
                { "Condition":[["A"],["B"],["C"]],
                  "Pr":[ 
                    [0,0,0,0,0,0,0.998,0,0,0,0,3.973E-4,0,1.329E-5,0.001,4.4E-203,0,0,0,0],
                    [0,0,0,0,0,0,0.954,0,0,0,0,0.003,0,0.042,2.66E-4,0,0,0,0,0],
                    [0,0,0.036,0,0,0,0.863,0,0,0,0,0.023,0,9.73E-4,0.034,0.036,0,0,0,0] ],
                  "Variable":{"Domain":{"Predef":"Protein"},"Name":"P25910__Pos184"},
                  "Nodetype":"CPT",
                  "Index":[0,1,2] },
                { "Condition":[["A"],["B"],["C"]],
                  "Pr":[
                    [0,0,0.081,0.185,0,0,0.733,0,0,0,0,3.83E-33,0,0,2.91E-4,1.65E-32,0,0,0,0],
                    [0,0,9.36E-4,0.001,0,0,0.977,0,0,0,0,1.59E-44,0,0,0.021,4.53E-44,0,0,0,0],
                    [0,0,3.61E-4,0.007,0,0,0.884,0,0,0,0,0.036,0,0.036,7.91E-4,0.036,0,0,0,0] ],
                  "Variable":{"Domain":{"Predef":"Protein"},"Name":"P25910__Pos186"},
                  "Nodetype":"CPT",
                  "Index":[0,1,2] },
                { "Condition":[["A"],["B"],["C"]],
                  "Pr":[
                    [0,0,0,0.167,0,0,0.766,0,0.007,0.008,0,0,0,4.45E-203,2.59E-7,0.0527,0,0,0,0],
                    [0,0,0,8.17E-6,0,0,0.996,0,2.29E-4,2.70E-4,0,0,0,0,1.64E-8,0.003,0,0,0,0],
                    [0.034,0,0,8.55E-5,0,0,0.721,0,0.028,0.026,0,0,0,0.038,0.113,0.036,0,0,0,0] ],
                  "Variable":{"Domain":{"Predef":"Protein"},"Name":"P25910__Pos306"},
                  "Nodetype":"CPT",
                  "Index":[0,1,2] } ],
              "Name":"P25910" },
          "Queries":[0,1,2,"Q704V1"] } }

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

### `<sequence>`

Examples of sequences as JSON:

    {"Seq":["T","N","G","A","T","C","C","N","N","G"],"Name":"S_1","Datatype":{"Predef":"DNA with N"}}
    {"Seq":["R","T","N","R","M","A","-","R","C","E","V","N","D","T","H","Y","L","-","R","V","D","R","N","V","E","R","M"],"Name":"S_2","Datatype":{"Predef":"Protein with gap"}}
    {"Seq":["U","G","C","U","G","U","G","G","U","A","U","C","G","G","U","U","A","G","A","A","U","C","G","G","G","A","A","C","C","C","A","G","A","G","C"],"Name":"S_3","Datatype":{"Predef":"RNA"}}
    {"Seq":["C","G","T","T","T","G","T","G","G","A","A","C","A","A","T","A","A","C","G","C","T","T","G","G","A","A","T","C","T","T"],"Name":"S_4","Datatype":{"Predef":"DNA"}}
    {"Seq":["W","L","M","F","N","Q","Y","T","T","S","C","H","A","A","D","G","M","I","K","F","N","Q","W","L","G","P","V","F","M","X","V"],"Name":"S_5","Datatype":{"Predef":"Protein with X"}}
    {"Seq":[false,true,true,true,false,true,true,true,true,true,true,false,true,false,true,false,false,false,false,true,true,true,false,false,false,true,true,true,true,true,false,true,false,false,false,true,false,true,false,false,true,false,true,true],"Name":"S_6","Datatype":{"Predef":"Boolean"}}
    {"Seq":["N"],"Name":"S_7","Datatype":{"Predef":"RNA with N"}}
    {"Seq":["L","P","A","Q","P","S","Q","N","Q","D","A","R","R"],"Name":"S_8","Datatype":{"Predef":"Protein"}}
    {"Seq":["G","T","A",null,"T","N",null,"T","G","N","A","A","A",null,"A","N",null,null,"G","A",null,"T",null,"C",null,null,null,"T"],"Name":"GS_1","Datatype":{"Predef":"DNA with N"}}
    {"Seq":[null,null,null,"S",null,null,null,null,"Q","E","L","S","T","W","A","-",null,null,"G","N","I","L"],"Name":"GS_2","Datatype":{"Predef":"Protein with gap"}}
    {"Seq":[],"Name":"GS_3","Datatype":{"Predef":"RNA"}}
    {"Seq":["A",null,"A","G",null,"T","A","C"],"Name":"GS_4","Datatype":{"Predef":"DNA"}}
    {"Seq":["F","V",null,null,"H","E","W","N","V","D","P",null,"M","C",null,null,"I","M","H","X","E","Q","E"],"Name":"GS_5","Datatype":{"Predef":"Protein with X"}}
    {"Seq":[false,null,null,false,false,null,null,false,true,true,false,true,true,true,true,false,null,false],"Name":"GS_6","Datatype":{"Predef":"Boolean"}}
    {"Seq":["U","C","N",null,"U","U","A","G","U","A","G",null,null,"C","A","G","N","A",null,null,"U","N","A","C","C",null,null,null,"C","G",null,"U"],"Name":"GS_7","Datatype":{"Predef":"RNA with N"}}
    {"Seq":["W",null,"S"],"Name":"GS_8","Datatype":{"Predef":"Protein"}}

### `<alignment>` (`dat.EnumSeq.Alignment`)

Example of alignment:

    { "Sequences":[
        {"Seq":[null,"G",null,"A","G",null,"G","G","A","G"],"Name":"AS_1"},
        {"Seq":[null,"G",null,"A","G",null,"G","G","T","C"],"Name":"AS_2"},
        {"Seq":[null,"G",null,"G","A",null,"G","T","A","A"],"Name":"AS_3"},
        {"Seq":["G","G","A","C","C",null,"A","T","C","G"],"Name":"AS_4"},
        {"Seq":["A","C",null,"T","G",null,"T","C","A","G"],"Name":"AS_5"}, 
        ...
      "Datatype":{"Predef":"DNA"} }


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

Currently, it is possible to represent either samples by "headers", or by "items" each of which has "features". The latter implies that data are indexed, not only by sample and header, but by sample, item and feature, forming a 3D tensor. 

    { "Headers":[<header1>,<header2>, ...],
      "Data":[ [ ... ] ] }

or

    { "Items":[<item1>,<item2>, ...],
      "Features":[<feature1>,<feature2>, ...],
      "Data":[ [ [ ... ] ] ] }

Note: the data matrix is first indexed by sample, then either header, or by item then feature.

#### Example

The example has 39 protein names as `Headers`, 
which are listed in the same order as the list of observations for each. 
In the example below, there is two observations for each protein name; 
the value `null` signifies absence of observation. 

    { "Headers":["A5ILB0","P08144","P29957","H2N0D4"],
      "Data":[[8.5,7.35,7.35,7],[9,9,null,8.5]]}

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











