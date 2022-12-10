/*
    bnkit -- software for building and using Bayesian networks
    Copyright (C) 2014  M. Boden et al.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package bn;

/**
 * Interface for BNodes that share parameters collected during training, i.e. counts.
 * Use by constructing BNode instances then "couple" one instance to another,
 * where the latter is the source for counts of the former.
 * Note that implies that "counts" are kept and accumulated only in the "source" node.
 *
 * @author mikael
 */
public interface TiedNode<T extends BNode> extends BNode {

    boolean tieTo(T master);

    /**
     * Return the TiedNode that is the "master" or "source" of this BNode.
     * @return master node (where parameters are sourced from learning), or null if this node is the master
     */
    T getMaster();
}
