package bn.ctmc;
import bn.ctmc.matrix.JCGap;
import bn.ctmc.matrix.JTTGap;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.assertEquals;

public class GapSubstModelTest {



    @Test
    void testConditionalProbabilitiesJCGap() {

        SubstModel JC = SubstModel.createModel("JC"); // only use to get domain
        double[][] IRM = {
                { 0, 0.333333, 0.333333, 0.333333},
                { 0, 0,        0.333333, 0.333333},
                { 0, 0,        0,        0.333333},
                { 0, 0,        0,        0       }
        };

        double mu = 0.01;
        double lambda = 0.02;
        // unnormalised JC in standard form i.e. mean rate is not 1 because we want to test the exact formula
        GapSubstModel JCGap = new JCGap(IRM, mu, lambda, true, true);

        double t = 0.05;
        double alpha = 0.3333333333;

        for (Object from : JC.getDomain().getValues()) { // Only go through non-gap states
            for (Object to : JC.getDomain().getValues()) {
                double probJCExact = transitionProbJCGapExactForm(from, to, mu, lambda, t, alpha);
                double probNoIndel = JCGap.getProb(from, to, t);
                assertEquals(probJCExact, probNoIndel, 1e-5, String.format("Conditional probabilities differ for %s -> %s", from, to));
            }
        }
    }

    /**
     * Calculation for extended JC69 model (equation 59 - Sup A1. Eddy and Rivas, 2008)
     */
    public double transitionProbJCGapExactForm(Object X, Object Y, double mu, double lambda, double t, double alpha) {

        double kronecker_delta = (X == Y) ? 1 : 0;

        return 0.25 * (lambda / (lambda + mu)) + (kronecker_delta - 0.25) * Math.exp(-(4 * alpha + mu) * t) + 0.25 * (mu / (lambda + mu)) * Math.exp(-(lambda + mu) * t);

    }

    @Test
    void testConditionalProbabilitiesJCNonGap() {

        SubstModel JC = SubstModel.createModel("JC"); // only use to get domain
        double[][] IRM = {
                { 0, 0.333333, 0.333333, 0.333333},
                { 0, 0,        0.333333, 0.333333},
                { 0, 0,        0,        0.333333},
                { 0, 0,        0,        0       }
        };
        double mu = 0;
        double lambda = 0;
        // un-normalised JC in standard form i.e. mean rate is not 1 because we want to test the exact formula
        GapSubstModel JCGap = new JCGap(IRM, mu, lambda, true, true);

        double t = 0.05;
        double alpha = 0.3333333333;

        for (Object from : JC.getDomain().getValues()) { // Only go through non-gap states
            for (Object to : JC.getDomain().getValues()) {
                double probJCExact = transitionProbJCExactForm(from, to, t, alpha);
                double probNoIndel = JCGap.getProb(from, to, t);
                assertEquals(probJCExact, probNoIndel, 1e-5, String.format("Conditional probabilities differ for %s -> %s", from, to));
            }
        }
    }

    public double transitionProbJCExactForm(Object X, Object Y, double t, double alpha) {
        double kronecker_delta = (X == Y) ? 1 : 0;

        return 0.25 + (kronecker_delta - 0.25) * Math.exp(-4 * alpha * t);
    }


    @Test
    void testConditionalProbabilitiesJCNoIndels() {

        // in this case we're comparing the standard normalised JC model to the JCGap with mu=lambda=0
        SubstModel JC = SubstModel.createModel("JC");
        GapSubstModel JC_no_indel = new JCGap(0, 0);
        double t = 0.5;

        for (Object from : JC.getDomain().getValues()) {
            for (Object to : JC.getDomain().getValues()) {
                double probJC = JC.getProb(from, to, t);
                double probNoIndel = JC_no_indel.getProb(from, to, t);
                assertEquals(probJC, probNoIndel, 1e-5, String.format("Conditional probabilities differ for %s -> %s", from, to));
            }
        }
    }


    @Test
    void testStationaryProbabilitiesJCNoIndels() {
        SubstModel JC = SubstModel.createModel("JC");
        GapSubstModel JC_no_indel = new JCGap(0, 0);

        for (Object state : JC.getDomain().getValues()) {
            double probJC = JC.getProb(state);
            double probNoIndel = JC_no_indel.getProb(state);
            assertEquals(probJC, probNoIndel, 1e-5, "Stationary probabilities differ for state: " + state);
        }
    }



    @Test
    void testConditionalProbabilitiesJTTNoIndels() {
        SubstModel JTT = SubstModel.createModel("JTT");
        GapSubstModel JTT_no_indel = new JTTGap(0, 0);
        double t = 0.5;

        for (Object from : JTT.getDomain().getValues()) {
            for (Object to : JTT.getDomain().getValues()) {
                double probJTT = JTT.getProb(from, to, t);
                double probNoIndel = JTT_no_indel.getProb(from, to, t);
                assertEquals(probJTT, probNoIndel, 1e-5, String.format("Conditional probabilities differ for %s -> %s", from, to));
            }
        }
    }

    @Test
    void testStationaryProbabilitiesJTTNoIndels() {
        SubstModel JTT = SubstModel.createModel("JTT");
        GapSubstModel JTT_no_indel = new JTTGap(0, 0);

        for (Object state : JTT.getDomain().getValues()) {
            double probJTT = JTT.getProb(state);
            double probNoIndel = JTT_no_indel.getProb(state);
            assertEquals(probJTT, probNoIndel, 1e-5, "Stationary probabilities differ for state: " + state);
        }
    }

    @Test
    void testGapProbabilities() {
        GapSubstModel JTT_no_indel = new JTTGap(0, 0);

        // Ensure that the gap probability is zero when there are no indels
        double gapProb = JTT_no_indel.getProb('-');
        assertEquals(0.0, gapProb, 1e-5, "Gap probability should be zero for JTT_no_indel");
    }

    @Test
    void testModelCopy() {
        GapSubstModel JTTGap = new JTTGap(0.05, 0.05);
        GapSubstModel JTTGapCopy = JTTGap.deepCopy();
        double t = 0.5;

        for (Object from : JTTGap.getDomain().getValues()) {
            for (Object to : JTTGap.getDomain().getValues()) {
                double probJTT = JTTGap.getProb(from, to, t);
                double probCopy = JTTGapCopy.getProb(from, to, t);
                assertEquals(probJTT, probCopy, 1e-5, String.format("Conditional probabilities differ for %s -> %s", from, to));
            }
        }
    }

}
