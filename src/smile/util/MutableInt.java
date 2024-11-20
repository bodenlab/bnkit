/*
 * Copyright (c) 2010-2021 Haifeng Li. All rights reserved.
 *
 * Smile is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Smile is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Smile.  If not, see <https://www.gnu.org/licenses/>.
 */

package smile.util;


/**
 * A mutable int wrapper. It is efficient as counter in HashMap.
 *
 * @author Haifeng Li
 */
public class MutableInt {

    /** The integer value. */
    public int value;

    /**
     * Constructor. The initial value is 1 since we're counting.
     */
    public MutableInt() {
        this(1);
    }

    /**
     * Constructor.
     *
     * @param value the initial value.
     */
    public MutableInt(int value) {
        this.value = value;
    }

    /**
     * Increment by one.
     * @return the update value.
     */
    public int increment() {
        return ++value;
    }

    /**
     * Increment.
     * @param x the operand.
     * @return the update value.
     */
    public int increment(int x) {
        return value += x;
    }

    /**
     * Decrement by one.
     * @return the update value.
     */
    public int decrement() {
        return --value;
    }

    /**
     * Decrement.
     * @param x the operand.
     * @return the update value.
     */
    public int decrement(int x) {
        return value -= x;
    }
}
