package smile.math.phylo;

import bn.ctmc.GapSubstModel;
import dat.EnumSeq;
import dat.Enumerable;
import dat.phylo.Tree;
import smile.math.Function;

/**
 * Function to compute the likelihood of an alignment given a tree and
 * a gap-augmented substitution model. Used for optimisation of indel rates.
 *
 * @author Sebastian
 */
public class AlnLikelihood implements Function {

    Tree tree;
    EnumSeq.Alignment<Enumerable> aln;
    double geometricSeqLenParam;
    Enumerable alpha;
    double[] F;
    double[][] IRM;

    public AlnLikelihood(Tree tree, EnumSeq.Alignment<Enumerable> aln, double geometricSeqLenParam,
                         Enumerable alpha, double[] F, double[][] IRM) {
        this.tree = tree;
        this.aln = aln;
        this.geometricSeqLenParam = geometricSeqLenParam;
        this.alpha = alpha;
        this.F = F;
        this.IRM = IRM;
    }

    /**
     * Get the log likelihood of a particular alignment given the tree and gap augmented substitution model.
     * This is used by a minimisation routine to find optimal params for the substitution model.
     * @param muLambda Assumes mu (deletion rate) and lambda (insertion rate) are equal
     * @return likelihood of the alignment given the tree + mu + lambda.
     */
    @Override
    public double f(double muLambda) {

        GapSubstModel model = new GapSubstModel(this.F, this.IRM, this.alpha, muLambda, muLambda);

        // trying to maximise the log likelihood
        return -1 * tree.calcAlnLikelihood(model, aln, geometricSeqLenParam, alpha);
    }
}
