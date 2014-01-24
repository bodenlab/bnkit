/**
 *
 */
package bn;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Enumerable variable class.
 *
 * @author mikael
 */
public class EnumVariable extends Variable<Enumerable> {

    public EnumVariable(Enumerable domain) {
        super(domain);
    }

    public EnumVariable(Enumerable domain, String name) {
        super(domain, name);
    }

    public int size() {
        return super.getDomain().order;
    }

    public int getIndex(Object value) {
        return super.getDomain().getIndex(value);
    }

    public static List<EnumVariable> toList(EnumVariable[] variables) {
        if (variables == null) {
            return null;
        }
        List<EnumVariable> list = new ArrayList<EnumVariable>(variables.length);
        for (EnumVariable var : variables) {
            list.add(var);
        }
        return list;
    }

}
