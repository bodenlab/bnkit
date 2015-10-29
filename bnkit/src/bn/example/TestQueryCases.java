package bn.example;

import bn.BNet;
import bn.BNode;
import bn.Distrib;
import bn.alg.CGTable;
import bn.alg.CGVarElim;
import bn.alg.Query;
import bn.file.BNBuf;
import bn.file.DataBuf;
import bn.prob.EnumDistrib;
import dat.EnumVariable;
import dat.Variable;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.*;

/**
 * Created by aesseb on 27-Oct-15.
 *
 * A class containing a variety of tests that can be performed on a trained Bayesian network.
 *
 * The user is provided five different tests that can be selected using setting:
 * Infer - given the network, data, query and sample number; perform inference (marginal distribution) and return a result for each row of data.
 *         Sampling is used to return the most accurate results based on the value set as max
 * MPE - given the network, data and query; using the joint distribution return the most likely state of query nodes given evidence
 * LLH_WI - given the network, data and query; perform inference and set each query node according to the result then calculate
 *          the likelihood of the model given the current state of the network
 * LLH_UD - given the network and data; set all nodes that have evidence then calculate the likelihood of the model given
 *          the current state of the network
 * Sample - given the network, data, query and sample number; perform inference and sample from the resulting distribution
 *          providing a set of results that can be plotted as a distribution. Effective for both real and discrete nodes
 *
 * @author aesseb
 */
public class TestQueryCases {

    private String setting = "";
    private int max = 0;
    private int states = 0; //how many rows of data are we testing?
    private Variable[] qVars; //how many variables in the query?

    /**
     *
     * @param bn_file - trained network in XML format
     * @param data_file - tab separated file with node names as header - not all nodes must have data
     * @param query - ; separated list of query node names
     * @param setting - choice of Infer, MPE, LLH_WI, LLH_UD, Sample
     * @param max - if using Sample, how many samples to generate
     */
    public TestQueryCases(String bn_file, String data_file, String query, String setting, int max) {

        this.setting = setting;
        this.max = max;

        BNet bn = BNBuf.load(bn_file); //load the network
        List<BNode> nodes = bn.getOrdered();
        Object[][] values = DataBuf.load(data_file, nodes); //load the data
        states = values.length;

        String[] queryS = query.split(";");
        Variable[] qVars = new Variable[queryS.length];
        for (int b = 0; b < queryS.length; b++) { //create a variable array of queries
            qVars[b] = (bn.getNode(queryS[b]).getVariable());
        }
        this.qVars = qVars;

        Variable[] vars = new Variable[nodes.size()];

        for (int k = 0; k < nodes.size(); k++) {
            vars[k] = nodes.get(k).getVariable();
        }

        testQuery(values, vars, bn, qVars); //based on original inputs, generate requested output

    }

    /**
     * Uses inputs to initialise network and create a Variable Elimination instance. Based on setting -
     * performs necessary tests, generates results and saves results
     *
     * @param values - data formatted as 2d array
     * @param vars - array of variables for ALL nodes in network
     * @param bn - loaded Bayesian network
     * @param qVars - array of query variables
     */
    public void testQuery(Object[][] values, Variable[] vars, BNet bn, Variable[] qVars) {

        Map<String, Map<Variable, Object>> store = new HashMap<>(); //to store results for Infer and MPE
        Map<String, Map<Variable, Object[]>> storeSample = new HashMap<>(); //sample case will return a different set of results
        Map<String, Double> storeLLH = new HashMap<>(); //LLH_WI and LLH_UD will return a different set of results
        Set<Variable> setNodes = new LinkedHashSet<>(); //FIXME this could change with null values?
        for (int i = 0; i < values.length; i++) {
            // set variables and keys according to observations
            String obs = ""; //record the observation for this row - used as key for storing results
            for (int j = 0; j < vars.length; j++) {
                BNode instantiate_me = bn.getNode(vars[j]);
                if (instantiate_me == null) {
                    System.out.println("Instantiate_me " + instantiate_me.getName() + " == null");
                }
                // check so that the observation is not null and instantiate_me is NOT a query variable
                if (values[i][j] != null && !Arrays.asList(qVars).contains(instantiate_me.getVariable())) {
                    // the node is instantiated to the value in the data set
                    instantiate_me.setInstance(values[i][j]);
                    setNodes.add(instantiate_me.getVariable()); //Record which nodes are instantiated
                    obs = obs + values[i][j] + ";"; //Record the current instantiation
                } else { // the observation is null or instantiate_me is a query variable and should remain uninstantiated
                    // the node is reset, i.e. un-instantiated
                    instantiate_me.resetInstance();
                }
            }

            CGVarElim ve = new CGVarElim();
            ve.instantiate(bn);

            Map<Variable, Object> result;
            switch (setting) {
                case "Infer":
                    Query q = ve.makeQuery(qVars);
                    CGTable cg = (CGTable) ve.infer(q);
                    result = infer(qVars, cg);
                    store.put(obs, result);
                    continue;
                case "MPE":
                    Query mpe = ve.makeMPE(qVars);
                    CGTable mpeOut = (CGTable) ve.infer(mpe);
                    result = MPE(qVars, mpeOut);
                    store.put(obs, result);
                    continue;
                case "LLH_WI":
                    Query ql = ve.makeQuery(qVars);
                    CGTable cgl = (CGTable) ve.infer(ql);
                    result = infer(qVars, cgl);
                    for (Map.Entry<Variable, Object> e : result.entrySet()) {
                        BNode instantiate_me = bn.getNode(e.getKey());
                        instantiate_me.setInstance(e.getValue());
                        obs = obs + instantiate_me.getVariable() + ";";
                        setNodes.add(instantiate_me.getVariable());
                    }
                    CGVarElim veLI = new CGVarElim();
                    veLI.instantiate(bn);
                    storeLLH.put(obs, veLI.likelihood());
                    continue;
                case "LLH_UD":
                    CGVarElim veL = new CGVarElim();
                    veL.instantiate(bn);
                    storeLLH.put(obs, veL.likelihood());
                    continue;
                case "Sample":
                    Query qs = ve.makeQuery(qVars);
                    CGTable cgs = (CGTable) ve.infer(qs);
                    Map<Variable, Object[]> res = sample(qVars, cgs);
                    storeSample.put(obs, res);
                    continue;
            }
        }
        switch(setting) {
            case "Infer":
                savePredictions(store, setNodes);
                break;
            case "MPE":
                savePredictions(store, setNodes);
                break;
            case "Sample":
                saveSample(storeSample, setNodes);
                break;
            case "LLH_WI":
                saveLLHs(storeLLH, setNodes);
                break;
            case "LLH_UD":
                saveLLHs(storeLLH, setNodes);
                break;
        }
    }

    /**
     * Calculate MPE for current set of query values
     *
     * @param qVars
     * @param mpeOut
     * @return a Map containing query variable as key and MP state as value
     */
    public Map<Variable, Object> MPE(Variable[] qVars, CGTable mpeOut) {

        Map<Variable, Object> results = new HashMap<>(qVars.length); //to store query variable results
        Variable.Assignment[] out = mpeOut.getMPE();
        for (Variable.Assignment a : out) {

            double sumPredic = 0.0;
            try { //If current query is discrete record the MPE
                EnumVariable v = (EnumVariable)a.var;
                results.put(a.var, a.val);
            } catch (ClassCastException e) { //If current query is real, sample from the distribution to to find most probable result
                //done by binning samples and reporting bin with highest count
                int b = 0;
                int max = 1000;
                while (b < max) {
                    Distrib d = (Distrib)a.val;
                    int nElem = 500; //number of elements
                    double minE = 0.0; //smallest value
                    double maxE = 5.0; //largest value FIXME - shouldn't be hard coded
                    double[] hist = new double[nElem];
                    double binSize = (maxE - minE)/nElem;
                    int x = 0;
                    double prediction;
                    while (x < 50000) { //sample from the distribution
                        double s = (Double) d.sample();
                        int bin = (int) ((s - minE) / binSize);
                        //Keep sampling from the same distribution until you get a reasonable result
                        if (bin < 0 ) { continue;}
                        else if (bin >= nElem) {continue;}
                        else {
                            hist[bin] += 1;
                            x++;
                        }
                    }
                    //**************************************************
                    //find the index for the bin with the most counts
                    int maxindex = 0;
                    for (int z = 0; z < hist.length; z++) {
                        double newnumber = hist[z];
                        if (newnumber > hist[maxindex])
                            maxindex = z;
                        //**************************************************
                    }
                    prediction = (double)(maxindex * binSize) + minE; //The bin with the highest counts
                    sumPredic += prediction;
                    b++;
                }
                double avg = sumPredic / b;
                results.put(a.var, avg);
            }
        }
        return results;
    }


    /**
     * Infer results for current set of queries
     *
     * @param qVars
     * @param cg
     * @return a Map containing query variable as key and MP state as value
     */
    public Map<Variable, Object> infer(Variable[] qVars, CGTable cg) {

        Map<Variable, Object> results = new HashMap<>(qVars.length); //to store query variable results
        for (Variable query : qVars) {
            double sumPredic = 0.0;
            try { //If query is discrete, sample max times and report the most common result
                EnumDistrib d = (EnumDistrib)cg.query(query);
                Map<Object, Object> counts = new HashMap<>();
                double a = 0;
                while (a < max) {
                    Object s = d.sample();
                    if (counts.containsKey(s)) {
                        counts.put(s, (Integer)counts.get(s) + 1);
                    } else {
                        counts.put(s, 1);
                    }
                    a++;
                }
                //Find most common result
                Integer mCount = 0;
                Object mKey = null;
                for (Map.Entry<Object, Object> c : counts.entrySet()) {
                    Object key = c.getKey();
                    Integer val = (Integer)c.getValue();
                    if (val > mCount) {
                        mCount = val;
                        mKey = key;
                    }
                }
                results.put(query, mKey);
            } catch (ClassCastException e) { //If current query is real, sample from the distribution to to find most probable result
                //done by binning samples and reporting bin with highest count
                int b = 0;
                while (b < max) {
                    Distrib d = cg.query(query);
                    double sum = 0;
                    int nElem = 500; //number of elements
                    double minE = 0.0; //smallest value
                    double maxE = 6.0; //largest value FIXME shouldn't be hard coded!
                    double[] hist = new double[nElem];
                    double binSize = (maxE - minE)/nElem;
                    int x = 0;
                    double prediction;
                    while (x < 10000) { //sample from the distribution
                        double s = (Double) d.sample();
                        int bin = (int) ((s - minE) / binSize);
                        //Keep sampling from the same distribution until you get a reasonable result
                        if (bin < 0 ) { continue;}
                        else if (bin >= nElem) {continue;}
                        else {
                            hist[bin] += 1;
                            x++;
                        }
                    }
                    //**************************************************
                    //find the index for the bin with the most counts
                    int maxindex = 0;
                    for (int z = 0; z < hist.length; z++) {
                        double newnumber = hist[z];
                        if (newnumber > hist[maxindex])
                            maxindex = z;
                        //**************************************************
                    }
                    prediction = (double)(maxindex * binSize) + minE; //The bin with the highest counts
                    sumPredic += prediction;
                    b++;
                }
                double avg = sumPredic / b;
                results.put(query, avg);
            }
        }
        return results;
    }

    /**
     * Similar to inference but all results are reported rather than the most common
     *
     * @param qVars
     * @param cg
     * @return a Map containing query variable as key and an array of sampled states as value
     */
    public Map<Variable, Object[]> sample(Variable[] qVars, CGTable cg) {

        Map<Variable, Object[]> results = new LinkedHashMap<>(qVars.length);
        for (Variable query : qVars) {
            int counter = 0;
            double sumPredic = 0.0;
            try { //If query is discrete
                EnumDistrib d = (EnumDistrib)cg.query(query);
                Object[] samples = new Object[max];
                counter = 0;
                while (counter < max) {
                    Object s = d.sample();
                    samples[counter] = s;
                    counter++;
                }
                results.put(query, samples);
            } catch (ClassCastException e) { //query is real
                counter = 0;
                Object[] samples = new Object[max];
                while (counter < max) {
                    Distrib d = cg.query(query);
                    double s = (Double) d.sample();
                    samples[counter] = s;
                    counter++;
                }
                results.put(query, samples);
            }
        }
        return results;
    }

    /**
     * Method to save results stored during sampling
     *
     * @param result
     * @param setNodes
     */
    public void saveSample(Map<String, Map<Variable, Object[]>> result, Set<Variable> setNodes) {

        int columns = setNodes.size() + qVars.length;
        Object[][] toPrint = new Object[(max*states)+1][columns];

        //Initialize the header of the output
        int s = 0;
        for (Variable vsn : setNodes) { //Record the nodes with information
            toPrint[0][s] = vsn;
            s++;
        }
        Object[] resKeys = new Object[qVars.length];
        int o = 0;
        for (Variable vres : qVars) { //Record the nodes that were queried
            toPrint[0][s] = vres;
            resKeys[o] = vres;
            s++;
            o++;
        }

        //Record the information for every set instance and the result from the query
        //In this case, formatted as a data frame that can be passed easily to R ggplot2
        int position = 1; //Need to track where in the bigger array we should start adding new info
        for (Map.Entry<String, Map<Variable, Object[]>> resEntry : result.entrySet()) {
            String key = resEntry.getKey();
            String[] obs = key.split(";");
            Map<Variable, Object[]> res = resEntry.getValue();
            for (int y = setNodes.size(); y < columns; y++) {
                Object[] curCol = res.get(resKeys[y - setNodes.size()]);
                for (int z = position; z < max + position; z++) {
                    for (int e = 0; e < setNodes.size(); e++) {
                        toPrint[z][e] = obs[e];
                    }
                    toPrint[z][y] = curCol[z-position]; //Record each result from the sampled array
                }
            }
            position += max;
        }

        try {
            PrintWriter writer = new PrintWriter("micro_data/" + setting + "_results.txt", "UTF-8");
            for (int j = 0; j < (max * states) + 1 ; j++) {
                for (int k = 0; k < columns; k++) {
                    if (k == columns - 1)
                        writer.write(toPrint[j][k] + "\n");
                    else
                        writer.write(toPrint[j][k] + "\t");
                }
            }
            writer.close();
        } catch (FileNotFoundException fnf) {
            System.out.println(fnf.getStackTrace());
        } catch (UnsupportedEncodingException use) {
            System.out.println(use.getStackTrace());
        }
    }

    /**
     * Method to save results stored during any setting except sampling
     *
     * @param result
     * @param setNodes
     */
    public void savePredictions(Map<String, Map<Variable, Object>> result, Set<Variable> setNodes) {

        int columns = setNodes.size() + qVars.length;
        Object[][] toPrint = new Object[states+1][columns];

        //Initialize the header of the output
        int s = 0;
        for (Variable vsn : setNodes) { //Record the nodes with information
            toPrint[0][s] = vsn;
            s++;
        }
        Object[] resKeys = new Object[qVars.length];
        int o = 0;
        for (Variable vres : qVars) { //Record the nodes that were queried
            toPrint[0][s] = vres;
            resKeys[o] = vres;
            s++;
            o++;
        }

        //Record the information for every set instance and the result from the query
        //In this case, formatted as a data frame that can be passed easily to R ggplot2
        int position = 1; //Need to track where in the bigger array we should start adding new info
        for (Map.Entry<String, Map<Variable, Object>> resEntry : result.entrySet()) {
            String key = resEntry.getKey();
            String[] obs = key.split(";");
            Map<Variable, Object> res = resEntry.getValue();

            for (int e = 0; e < setNodes.size(); e++) { //Record the set instances
                toPrint[position][e] = obs[e];
            }
            for (int col = setNodes.size(); col < columns; col++) { //Record the query results
                toPrint[position][col] = res.get(resKeys[col - setNodes.size()]);
            }
            position ++;
        }

        try {
            PrintWriter writer = new PrintWriter("micro_data/" + setting + "_results.txt", "UTF-8");
            for (int j = 0; j < states + 1 ; j++) {
                for (int k = 0; k < columns; k++) {
                    if (k == columns - 1)
                        writer.write(toPrint[j][k] + "\n");
                    else
                        writer.write(toPrint[j][k] + "\t");
                }
            }
            writer.close();
        } catch (FileNotFoundException fnf) {
            System.out.println(fnf.getStackTrace());
        } catch (UnsupportedEncodingException use) {
            System.out.println(use.getStackTrace());
        }
    }


    public void saveLLHs(Map<String, Double> result, Set<Variable> setNodes) {

        int columns = setNodes.size() + 1; //In this case only, all original and queried nodes become 'set'
        Object[][] toPrint = new Object[states+1][columns];

        //Initialize the header of the output
        int s = 0;
        for (Variable vsn : setNodes) { //Record the nodes with information
            toPrint[0][s] = vsn;
            s++;
        }
        toPrint[0][setNodes.size()] = "Likelihood";

        //Record the information for every set instance and the result from the query
        //In this case, formatted as a data frame that can be passed easily to R ggplot2
        int position = 1; //Need to track where in the bigger array we should start adding new info
        for (Map.Entry<String, Double> resEntry : result.entrySet()) {
            String key = resEntry.getKey();
            String[] obs = key.split(";");
            double value = resEntry.getValue();

            for (int e = 0; e < setNodes.size(); e++) { //Record the set instances
                toPrint[position][e] = obs[e];
            }

            toPrint[position][setNodes.size()] = value;

            position ++;
        }

        try {
            PrintWriter writer = new PrintWriter("micro_data/" + setting + "_results.txt", "UTF-8");
            for (int j = 0; j < states + 1 ; j++) {
                for (int k = 0; k < columns; k++) {
                    if (k == columns - 1)
                        writer.write(toPrint[j][k] + "\n");
                    else
                        writer.write(toPrint[j][k] + "\t");
                }
            }
            writer.close();
        } catch (FileNotFoundException fnf) {
            System.out.println(fnf.getStackTrace());
        } catch (UnsupportedEncodingException use) {
            System.out.println(use.getStackTrace());
        }
    }

}
