package smile.math.phylo;

import bn.ctmc.GapSubstModel;
import dat.EnumSeq;
import dat.Enumerable;
import dat.phylo.Tree;
import smile.math.Function;

public class AlnLikelihood implements Function {

    Tree tree;
    EnumSeq.Alignment<Enumerable> aln;
    double geometricSeqLenParam;
    Enumerable alpha;
    GapSubstModel model;

    public AlnLikelihood(Tree tree, EnumSeq.Alignment<Enumerable> aln, double geometricSeqLenParam,
                         Enumerable alpha, GapSubstModel model) {
        this.tree = tree;
        this.aln = aln;
        this.geometricSeqLenParam = geometricSeqLenParam;
        this.alpha = alpha;
        this.model = model;
    }

    /**
     * Get the log likelihood of a particular alignment given the tree and gap augmented substitution model.
     * This is used by a minimisation routine to find optimal params for the substitution model.
     * @param mu_lambda Assumes mu (deletion rate) and lambda (insertion rate) are equal
     * @return likelihood of the alignment given the tree + mu + lambda.
     */
    @Override
    public double f(double mu_lambda) {

        GapSubstModel model = new GapSubstModel(this.model.getOrigF(), this.model.getIRM(),
                this.model.getOrigAlpha(), mu_lambda, mu_lambda);

        // trying to maximise the log likelihood
        return -1 * tree.calcAlnLikelihood(model, aln, geometricSeqLenParam, alpha);
    }
}
