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

package dat;

import json.JSONObject;

import java.util.HashMap;
import java.util.Map;

/**
 * Continuous domain definition. 
 * For checking validity of values for variables that belong to this domain.
 * @author mikael
 */
public class Continuous implements Domain {

    public boolean isValid(Object value) {
        try {
            Double x = (Double) value;
            return true;
        } catch (ClassCastException e) {
            return false;
        }
    }

    /**
     * Two Continuous domains are equal by definition
     * @param other
     * @return
     */
    @Override
    public boolean equals(Object other) {
        try {
            Continuous x = (Continuous) other;
            return true;
        } catch (ClassCastException e) {
            return false;
        }
    }

    @Override
    public JSONObject toJSON() {
        JSONObject jobj = new JSONObject();
        jobj.put("Predef", "Real");
        return jobj;
    }

    private static Map<String, Continuous> predef = new HashMap<>();
    private static Map<Continuous, String> predef_reverse = new HashMap<>();

    private static Continuous real = new Continuous();
    // Instantiating the static map
    static
    {
        predef.put("Real", real);
        predef_reverse.put(real, "Real");
    }

    public static Continuous fromJSON(JSONObject json) {
        String defID = json.optString("Predef", null);
        if (defID != null) {
            Continuous match = predef.get(defID);
            if (match == null) {
                throw new RuntimeException("Invalid predefined domain name: \"" + defID + "\". Please check dat.Enumerable");
            }
            return match;
        }
        return real;
    }

}
