package util;

import java.util.Arrays;


public abstract class Binner {

    protected final int numBins;
    protected double[] binEdges;
    protected double[] sortedData;

    public Binner(int numBins) {
        this.numBins = numBins;
        binEdges = new double[numBins + 1];
    }

    public void fit(double[] data, boolean dereplicate) {
        sortedData = checkDataValid(data, dereplicate);
        binEdges = computeBinEdges(sortedData, numBins);
    }

    protected static double[] checkDataValid(double[] data, boolean dereplicate) {
        if (data == null || data.length == 0)
            throw new IllegalArgumentException("Data must not be empty");

        double[] sorted = data.clone();
        if (dereplicate) {
            sorted = Arrays.stream(sorted).distinct().toArray();
        }
        Arrays.sort(sorted);

        return sorted;
    }

    public void fit(double[] data) {
        if (data == null || data.length == 0)
            throw new IllegalArgumentException("Data must not be empty");

        double[] sorted = data.clone();
        Arrays.sort(sorted);
        sortedData = sorted;
        binEdges = computeBinEdges(sortedData, numBins);
    }

    protected abstract double[] computeBinEdges(double[] sorted, int nBuckets);

    public int transform(double value) {
        return bucket(value, binEdges);
    }

    public int[] transform(double[] data) {
        int[] result = new int[data.length];
        for (int i = 0; i < data.length; i++)
            result[i] = bucket(data[i], binEdges);
        return result;
    }

    protected static int bucket(double value, double[] binEdges) {
        int lo = 0, hi = binEdges.length - 2;  // hi = last valid bucket index
        while (lo < hi) {
            int mid = (lo + hi) / 2;
            if (value < binEdges[mid + 1]) hi = mid;
            else lo = mid + 1;
        }
        return lo;
    }

    public double[] getBinEdges() {
        return this.binEdges;
    }


    /**
     * A utility class for splitting data into even(ish) sized chunks.
     * Based off the numpy.array_split function. Here it is implemented in a way
     * that the object jointly stores the data, split indices and also thresholds
     * between each of the splits.
     */
    public static class EqualWidthBinner extends Binner {

        double[] groupMeans;
        int[] splitIndices;

        public EqualWidthBinner( int numSections) {
            super(numSections);
        }

        public void fit(double[] data) {
            if (data == null || data.length == 0) {
                throw new IllegalArgumentException("Data must not be empty");
            }

            double[] sorted = data.clone();
            Arrays.sort(sorted);
            sortedData = sorted;
            this.splitIndices = calcSplitIndices(sortedData, numBins);
            binEdges = computeBinEdges(sortedData, numBins);
        }

        @Override
        protected double[] computeBinEdges(double[] sorted, int nBuckets) {

            double[] binEdges = new double[numBins + 1];
            double[] groupMeans = new double[numBins];
            for (int i = 0; i < numBins; ++i) {

                int start = splitIndices[i];
                int end = splitIndices[i + 1];

                double sum = 0.0;
                for (int j = start; j < end; ++j) {
                    sum += sortedData[j];
                }
                groupMeans[i] = sum / (end - start);
            }
            this.groupMeans = groupMeans;

            for (int i = 0; i < numBins - 1; ++i) {
                int end = splitIndices[i + 1];
                double currentGroupMax = sortedData[end - 1];
                double nextGroupMin = sortedData[end];
                binEdges[i + 1] = (currentGroupMax + nextGroupMin) / 2.0;
            }

            binEdges[0] = Double.NEGATIVE_INFINITY;
            binEdges[numBins] = Double.POSITIVE_INFINITY;

            return binEdges;

        }

        /**
         * Based off numpy.array_split. For an array of length l that should be
         * split into n sections, it returns an array which defines the start of each
         * new section.
         */
        private int[] calcSplitIndices(double[] inputData, int numSections)  {

            int numData = inputData.length;
            int numEachSection = Math.floorDiv( numData, numSections);
            int extras = numData % numSections;

            int totalSections = 1 + extras + (numSections - extras);
            int [] sectionSizes = new int[totalSections];
            sectionSizes[0] = 0;

            // section_sizes = ([0] + extras * [n_each_section + 1] + (n_sections - extras) * [n_each_section])
            for (int j = 1; j < sectionSizes.length; j++) {
                if (j <= extras) {
                    sectionSizes[j] = numEachSection + 1;
                } else {
                    sectionSizes[j] = numEachSection;
                }
            }

            int[] divPoints = new int[totalSections];
            int cumulative = 0;
            for (int i = 0; i < sectionSizes.length; ++i) {
                cumulative += sectionSizes[i];
                divPoints[i] = cumulative;
            }

            return divPoints;
        }



        public static void main(String[] args) {
            double[] data = {5.0, 2.0, 9.0, 1.0, 3.0};
            int numSections = 4;
            Binner splitter = new EqualWidthBinner(numSections);
            splitter.fit(data);
            System.out.println("Data: " + java.util.Arrays.toString(splitter.sortedData));
            System.out.println("Thresholds: " + java.util.Arrays.toString(splitter.getBinEdges()));
            System.out.println("Bucket for 2.5: " + splitter.transform(2.49));
        }

    }

    public static class QuantileBinner extends Binner {

        private final boolean zeroInflated;

        public QuantileBinner(int numSections, boolean zeroInflated) {
            super(numSections);
            this.zeroInflated = zeroInflated;
        }

        public void fit(double[] data) {
            fit(data, false);
        }

        public void fit(double[] data, boolean dereplicate) {
            sortedData = checkDataValid(data, dereplicate);

            if (zeroInflated) {
                double[] nonZeroData = Arrays.stream(sortedData).filter(d -> d > 0).toArray();
                binEdges = computeBinEdges(nonZeroData, numBins - 1);
            } else {
                binEdges = computeBinEdges(sortedData, numBins);
            }
        }

        @Override
        protected double[] computeBinEdges(double[] sorted, int nBuckets)  {

            double[] binEdges = new double[nBuckets + 1];
            for (int i = 0; i <= nBuckets; i++) {
                double p = (double) i / nBuckets;  // 0.0 → 1.0
                binEdges[i] = percentile(sorted, p);
            }

            // Open-ended boundaries to handle out-of-range values at inference
            binEdges[0] = Double.NEGATIVE_INFINITY;
            binEdges[nBuckets] = Double.POSITIVE_INFINITY;

            return binEdges;
        }

        @Override
        public int transform(double value) {
            if (zeroInflated && value == 0.0) {
                return 0;  // Assign zero values to the first bucket
            }

            int bucketIndex = bucket(value, binEdges);
            if (zeroInflated) {
                return bucketIndex + 1;  // Shift by 1 to account for zero bucket
            } else {
                return bucketIndex;
            }
        }

        @Override
        public int[] transform(double[] data) {
            int[] result = new int[data.length];
            for (int i = 0; i < data.length; i++)
                result[i] = transform(data[i]);
            return result;
        }


        private static double percentile(double[] sorted, double p) {
            if (p == 0.0) return sorted[0];
            if (p == 1.0) return sorted[sorted.length - 1];

            double index = p * (sorted.length - 1);
            int lo = (int) index;
            double frac = index - lo;
            return sorted[lo] + (sorted[lo + 1] - sorted[lo]) * frac;  // i + (j-i)*fraction where i < j
        }

        public static void main(String[] args) {
            double[] data = {0, 0, 9.0, 1.0, 3.0};
            int numSections = 4;
            Binner splitter = new QuantileBinner(numSections, true);
            splitter.fit(data);
            System.out.println("Data: " + java.util.Arrays.toString(splitter.sortedData));
            System.out.println("Thresholds: " + java.util.Arrays.toString(splitter.getBinEdges()));
            System.out.println("Bucket for 0.1: " + splitter.transform(0.1));
        }

    }
}
