package dat.phylo;

/**
 * Intended to specify a module that will allow for a phylogenetic tree to be decorated
 * with whatever values, predictions, distributions etc each item to be a generic datatype E.
 * @param <E> what the tree is decorated with
 */
public interface TreeDecor<E> {
    E getDecoration(int idx);
    void decorate(TreeInstance ti);
}

