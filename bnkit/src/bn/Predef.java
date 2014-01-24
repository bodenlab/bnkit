package bn;

import java.util.ArrayList;
import java.util.List;

/**
 * This class defines a number of conventional variable and node types.
 * For a variable and node to be specified in a XML file to be loaded or saved, 
 * their types must be defined and mapped in this class. Furthermore, they need
 * to be recognised by methods that check if they accept parameters etc.
 * @author mikael
 */
public class Predef {

    public Predef() {
    }

    public static boolean isParameterised(String typename) {
        if (typename.equalsIgnoreCase("String") || (typename.equalsIgnoreCase("Number"))) {
            return true;
        } else {
            return false;
        }
    }

    public static String[] getVariableTypes() {
        return new String[]{"Boolean", "String", "Number", "Real", "Amino acid", "Nucleic acid"};
    }

    public static boolean isEnumerable(String typename) {
        if (typename.equalsIgnoreCase("Real")) {
            return false;
        } else {
            return true;
        }
    }

    public static String getBNodeType(String var_type) {
        if (isEnumerable(var_type)) {
            return "CPT";
        } else {
            return "GDT";
        }
    }

    @SuppressWarnings("rawtypes")
    public static Variable getVariable(String varname, String typename, String params) {
        if (typename.equalsIgnoreCase("Boolean")) {
            return Predef.Boolean(varname);
        }
        if (typename.equalsIgnoreCase("String")) {
            String[] values = params.split(";");
            return Predef.Nominal(values, varname);
        }
        if (typename.equalsIgnoreCase("Number")) {
            try {
                Integer max = Integer.parseInt(params);
                return Predef.Number(max, varname);
            } catch (NumberFormatException e) {
                e.printStackTrace();
                return null;
            }
        }
        if (typename.equalsIgnoreCase("Amino acid")) {
            return Predef.AminoAcid(varname);
        }
        if (typename.equalsIgnoreCase("Nucleic acid")) {
            return Predef.NucleicAcid(varname);
        }
        if (typename.equalsIgnoreCase("Real")) {
            return Predef.Real(varname);
        }
        throw new RuntimeException("Invalid specification of variable");
    }

    @SuppressWarnings("rawtypes")
    public static BNode getBNode(Variable var, List<Variable> parents, String type) {
        if (parents == null) {
            parents = new ArrayList<>();
        }
        try {
            if (type.equalsIgnoreCase("CPT")) {
                List<EnumVariable> elist = new ArrayList<>();
                for (Variable v : parents) {
                    elist.add((EnumVariable) v);
                }
                return new CPT((EnumVariable) var, elist);
            } else if (type.equalsIgnoreCase("GDT")) {
                List<EnumVariable> elist = new ArrayList<>();
                for (Variable v : parents) {
                    elist.add((EnumVariable) v);
                }
                return new GDT((Variable<Continuous>) var, elist);
            }
            return null;
        } catch (ClassCastException e) {
            return null;
        }
    }

    @SuppressWarnings("rawtypes")
    public static String getType(Variable var) {
        String type = var.getPredef();
        if (type != null) {
            return type;
        }
        return null;
    }

    public static String getType(BNode node) {
        return node.getType();
    }

    /**
     * Get the object representation of a string, for a specified variable.
     *
     * @param var variable (which needs to be a predefined type (@see bn.Variable#setPredef(String))
     * @param vstr string representation of value (as it appears in say a file)
     * @return the Object instance
     */
    public static Object getObject(Variable var, String vstr) {
        try {
            if (var.getPredef().equals("Boolean")) {
                return Boolean.parseBoolean(vstr);
            } else if (var.getPredef().equals("Nucleic acid")) {
                Character ch = new Character(vstr.charAt(0));
                if (var.getDomain().isValid(ch)) {
                    return ch;
                }
            } else if (var.getPredef().equals("Amino acid")) {
                Character ch = new Character(vstr.charAt(0));
                if (var.getDomain().isValid(ch)) {
                    return ch;
                }
            } else if (var.getPredef().equals("String")) {
                if (var.getDomain().isValid(vstr)) {
                    return vstr;
                }
            } else if (var.getPredef().equals("Number")) {
                Integer n = Integer.parseInt(vstr);
                if (var.getDomain().isValid(n)) {
                    return n;
                }
            } else if (var.getPredef().equals("Real")) {
                Double y = Double.parseDouble(vstr);
                if (var.getDomain().isValid(y)) {
                    return y;
                }
            }
        } catch (NumberFormatException e) {
            return null;
        }
        return null;
    }

    public static EnumVariable Boolean(String name) {
        EnumVariable var = new EnumVariable(Enumerable.bool, name);
        var.setPredef("Boolean");
        return var;
    }

    public static EnumVariable Boolean() {
        return Predef.Boolean("Bool");
    }

    public static EnumVariable NucleicAcid(String name) {
        EnumVariable var = new EnumVariable(Enumerable.nacid, name);
        var.setPredef("Nucleic acid");
        return var;
    }

    public static EnumVariable NucleicAcid() {
        return Predef.NucleicAcid("NT");
    }

    public static EnumVariable AminoAcid(String name) {
        EnumVariable var = new EnumVariable(Enumerable.aacid, name);
        var.setPredef("Amino acid");
        return var;
    }

    public static EnumVariable AminoAcid() {
        return Predef.AminoAcid("AA");
    }

    public static EnumVariable Nominal(String[] values, String name) {
        Enumerable domain = new Enumerable(values);
        EnumVariable var = new EnumVariable(domain, name);
        StringBuilder sbuf = new StringBuilder("");
        for (Object v : values) {
            sbuf.append(v).append(";");
        }
        var.setParams(sbuf.toString());
        var.setPredef("String");
        return var;
    }

    public static EnumVariable Nominal(String[] values) {
        return Predef.Nominal(values, "Nom");
    }

    public static EnumVariable Number(int max, String name) {
        Enumerable domain = new Enumerable(max);
        EnumVariable var = new EnumVariable(domain, name);
        var.setName("Number");
        var.setParams("" + max);
        return var;
    }

    public static EnumVariable Number(int max) {
        return Predef.Number(max, "Num");
    }

    public static Variable<Continuous> Real(String name) {
        Continuous domain = new Continuous();
        Variable<Continuous> var = new Variable<>(domain, name);
        var.setPredef("Real");
        return var;
    }

    public static Variable<Continuous> Real() {
        return Predef.Real("Rl");
    }

}
