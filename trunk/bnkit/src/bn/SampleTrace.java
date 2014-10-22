package bn;

import bn.factor.Factor;
import dat.EnumVariable;
import dat.Variable;
import dat.Variable.Assignment;
import bn.alg.ApproxInference;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;


/**
 * A structure used by Approximate inference to query a hybrid network. Allows
 * for counts and samples to be stored together and then processed into
 * appropriate structures based on the query.
 *
 * Records each 'state' of the Markov chain as it proceeds and in this way
 * records the influence of non-evidenced variables.
 *
 * @author Alex
 * @author Mikael (secondary)
 */
public class SampleTrace {

    private int[] enumRecord; // all observations of enumerable variables in query, encoded as key-indices
    private CountTable enumCounts; // the counts of enumerable variables in query
    private Map<BNode, SampleTable> nonEnumSamples;
    protected final int nQuery;
    protected final int nSamples;
    protected final List<Variable> qvars;       // All query variables
    private final List<EnumVariable> enumVars = new ArrayList<>();  // All enumerable query variables
    private final List<BNode> enumNodes = new ArrayList<>();  // All enumerable query nodes
    private final List<BNode> nonEnumNodes = new ArrayList<>();  // All non-enumerable query nodes
    private final List<Variable> nonEnumVars = new ArrayList<>();  // All non-enumerable query variables
    public Map<Variable, Map<Integer, List<Object>>> nonEnumTables = null;
    private boolean allNonEnumerable = false;
    private boolean allEnumerable = false;
    private int round = 0;
    private final Object[] key;
    
    /**
     * Create a container for recording samples of mixed types (enumerable and non-enumerable).
     * @param query the nodes which are recorded
     */
    public SampleTrace(List<BNode> query, int nSamples) {
        this.qvars = new ArrayList<>(query.size());     // all query variables
        this.nQuery = query.size();                     // the number of query variables
        this.nonEnumSamples = new HashMap<>();
        this.enumRecord = new int[nSamples];
        this.nSamples = nSamples;
        for (BNode node : query) {
            Variable var = node.getVariable();
            this.qvars.add(var);
            try {
                EnumVariable e = (EnumVariable) var;
                enumVars.add(e);
                enumNodes.add(node);
            } catch (ClassCastException e) {
                nonEnumVars.add(var);
                nonEnumNodes.add(node);
            }
        }
        
        if (enumNodes.isEmpty()) {
            allNonEnumerable = true;
            key = null;
        } else {
            key = new Object[enumNodes.size()];
            this.enumCounts = new CountTable(enumVars);
        }
        if (nonEnumNodes.isEmpty()) {
            allEnumerable = true;
        } else {
            for (BNode node : nonEnumNodes) 
                this.nonEnumSamples.put(node, new SampleTable(enumVars));                     // map of query samples for all non-enumerable variables CONDITIONED on enumerables
        }
    }

    /**
     * Add the current instances of the query variables to the appropriate list.
     * All nodes are assumed to be instantiated.
     */
    public void count() {
        int ecnt = 0;
        int key_index = 0;
        if (!allNonEnumerable) {
            for (BNode node : enumNodes)
                key[ecnt ++] = node.getInstance();
            key_index = enumCounts.getIndex(key);
            enumRecord[round] = key_index;
            enumCounts.count(key_index);
        }
        if (!allEnumerable) {
            if (round > 0) {
                for (Map.Entry<BNode, SampleTable> entry : nonEnumSamples.entrySet()) {
                    BNode node = entry.getKey();
                    SampleTable samples = entry.getValue();
                    if (allNonEnumerable) {
                    	samples.count(node.getInstance());
                    } else {
                    	samples.count(enumRecord[round - 1], node.getInstance());
                    }
                }
            }
        }
        round ++;
    }

    /**
     * Construct a factor from a CountTable for any number of enumerable variables, 
     * and a SampleTable for any number of non-enumerable variables CONDITIONED on the 
     * former enumerable variables.
     * @param counts
     * @param nonEnumSamples 
     */
    private Factor makeFactor(CountTable counts, Map<BNode, SampleTable> nonEnumSamples) {
        Factor f = new Factor(qvars);
        for (int index = 0; index < counts.table.getSize(); index ++) {
            double value = counts.get(index);
            f.setFactor(index, value);
            JDF jdf = new JDF(nonEnumVars);
            f.setJDF(index, jdf);
            for (Entry<BNode, SampleTable> sample : nonEnumSamples.entrySet()) {
                BNode node = sample.getKey();
                Variable var = node.getVariable();
                SampleTable table = sample.getValue();
                Distrib d = node.makeDistrib(table.getAll(index));
                jdf.setDistrib(d, var);
            }
        }
        return f;
    }
    
//    /**
//     * Construct a factor from a CountTable for any number of enumerable variables, 
//     * and a SampleTable for any number of non-enumerable variables CONDITIONED on the 
//     * former enumerable variables.
//     * @param counts
//     * @param nonEnumSamples 
//     */
//    private factor.AbstractFactor makeFactor(CountTable counts, Map<BNode, SampleTable> nonEnumSamples) {
//        Factor f = new Factor(qvars);
//        for (int index = 0; index < counts.table.getSize(); index ++) {
//            double value = counts.get(index);
//            f.setFactor(index, value);
//            JDF jdf = new JDF(nonEnumVars);
//            f.setJDF(index, jdf);
//            for (Entry<BNode, SampleTable> sample : nonEnumSamples.entrySet()) {
//                BNode node = sample.getKey();
//                Variable var = node.getVariable();
//                SampleTable table = sample.getValue();
//                Distrib d = node.makeDistrib(table.getAll(index));
//                jdf.setDistrib(d, var);
//            }
//        }
//        return f;
//    }
    
    /**
     * Construct a factor from a map of SampleTables (one for each non-enumerable variable) 
     * @param nonEnumSamples 
     */
    private Factor makeFactor(Map<BNode, SampleTable> nonEnumSamples) {
        Factor f = new Factor(qvars);
        f.setFactor(1.0);
        JDF jdf = new JDF(nonEnumVars);
        f.setJDF(jdf);
        for (Entry<BNode, SampleTable> entry : nonEnumSamples.entrySet()) {
            BNode node = entry.getKey();
            Variable var = node.getVariable();
            SampleTable table = entry.getValue();
            Distrib d = node.makeDistrib(table.getAll());
            jdf.setDistrib(d, var);
        }
        return f;
    }
    
    /**
     * Construct a factor with enumerable variables from a count table.
     * @param counts 
     */
    private Factor makeFactor(CountTable counts) {
        Factor f = new Factor(qvars);
        for (Entry<Integer, Double> entry : counts.table.getMapEntries()) {
            int index = entry.getKey();
            double value = entry.getValue();
            f.setFactor(index, value);
        }
        return f;
    }
    
    public Factor getFactor() {
        Factor f;
        if (allNonEnumerable) 
            f = makeFactor(this.nonEnumSamples);
        else if (allEnumerable)
            f = makeFactor(this.enumCounts);
        else
            f = makeFactor(this.enumCounts, this.nonEnumSamples);
        return f;
    }
}
