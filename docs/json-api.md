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

### Run reconstruction: `<command>` is `"Recon"`

`<params>` is 

`{"Tree":<tree>,"Alignment":<alignment>, <optional>}`

or 

`{"POGTree":<POGtree>,"Ancestors":[<POG1>,<POG2>, ...], <optional>}`

and

`<optional>` is

`"Inference":"Joint"` (default) or `"Inference":"Marginal","Ancestor":<anc-ID>`

`"Indels":<indel-method>` where `<indel-method>` is `"BEP"` (default), `"BEML"`, `"SICP"`, `"SICML"`, `"PSP"`, or `"PSML"`.

`"Model":<subst-model>` where `<subst-model>` is `"JTT"` (default), `"Dayhoff"`, `"LG"`, `"WAG"`, `"JC"`, or `"Yang"`.

`"Rates":[<rate-pos1>,<rate-pos2>, ...]`

Result is
    
    { "Prediction":{"Input":{<input-JSON>},
                    "Ancestors":[<list-of-POGs>] 
                    "Datatype":"Prediction"}}}

## Data structures

### `<tree>`

    { "Branchpoints":<number-of-idxs>,
      "Labels":[<list-of-labels>],
      "Parents":[<list-of-parent-idxs>],     // -1 represents no parent
      "Distances":[<list-of-dists>] }

#### Example

Newick definition

    ((A:0.6,((B:3.3,(C:1.0,D:2.5)cd:1.8)bcd:5,((E:3.9,F:4.5)ef:2.5,G:0.3)efg:7)X:3.2)Y:0.5,H:1.1)I:0.2

JSON definition (note internal nodes are re-numbered, left-to-right, depth-first)

    {"Parents":[-1,0,1,1,3,4,4,6,6,3,9,10,10,9,0],"Labels":["0","1","A","2","3","B","4","C","D","5","6","E","F","G","H"],"Distances":[0,0.5,0.6,3.2,5,3.3,1.8,1,2.5,7,2.5,3.9,4.5,0.3,1.1],"Branchpoints":15}

### `<alignment>`


### `<POG>`


### '<POGTree>'








