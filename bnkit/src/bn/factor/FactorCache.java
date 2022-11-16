package bn.factor;

import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

public class FactorCache {

    public final ConcurrentHashMap<Integer, CountedFactorMap> cache;

    public FactorCache() {
        cache = new ConcurrentHashMap();
    }

    public int size() {
        return cache.size();
    }

    public void reportCache() {
        int cnt = 1;
        for (Map.Entry<Integer, CountedFactorMap> entry : cache.entrySet()) {
            AbstractFactor.FactorMap fmap = entry.getValue().fmap;
            System.out.println((cnt ++) + "\t" + entry.getKey() + "\t" + entry.getValue().fmap.size() + "\t" + entry.getValue().accesses);
        }
    }

    public AbstractFactor.FactorMap get(int key) {
        CountedFactorMap cfmap = cache.get(key);
        if (cfmap != null)
            return cfmap.access();
        return null;
    }

    public void put(int key, AbstractFactor.FactorMap map) {
        CountedFactorMap cfmap = cache.get(key);
        if (cfmap == null)
            cache.put(key, new CountedFactorMap(map));
    }

    private class CountedFactorMap {
        AbstractFactor.FactorMap fmap;
        int accesses = 0;

        CountedFactorMap(AbstractFactor.FactorMap fmap) {
            this.fmap = fmap;
            this.accesses = 0;
        }
        AbstractFactor.FactorMap access() {
            accesses += 1;
            return fmap;
        }
    }
}
