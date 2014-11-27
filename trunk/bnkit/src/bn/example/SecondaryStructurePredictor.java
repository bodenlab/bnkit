package bn.example;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Iterator;

import dat.Domain;
import dat.EnumVariable;
import dat.Enumerable;
import bn.prob.EnumDistrib;
import bn.alg.CGTable;
import bn.alg.Query;
import bn.alg.VarElim;
import bn.node.CPT;
import bn.BNet;
import bn.Predef;
/**
 * used for secondary structure prediction
 * with algiment 
 * @author wangyufei
 */

public class SecondaryStructurePredictor {
	
	final static private int windowSize = 5;
	final static private String fileName = "data/points.csv";
	final static private double trainDataRadio = 0.6f;
	
	public static void main(String[] args) {
		EnumVariable classNode = Predef.SecondaryStructure("secondary structure class");
		BNet bn = new BNet();
	
		//train
		List<EnumVariable> features = new LinkedList<EnumVariable>();
		List<Character[]> data = loadData("data/points.csv",windowSize);
		Collections.shuffle(data);
		int cutOffIndex = (int)(data.size() * trainDataRadio);
		
		List<Character[]> training = data.subList(0, cutOffIndex);
		List<Character[]> testing = data.subList(cutOffIndex, data.size());
		
		for(int i = 0; i < windowSize; i++) {
			features.add(Predef.AminoAcid("feature" + String.valueOf(i)));
		}
		double[][][] trainResult = learn(training,features,classNode);
		
		CPT head = new CPT(classNode);
		EnumDistrib classDistrib = new EnumDistrib(classNode.getDomain());
		double[] classCount = new double[classNode.size()];
		
		for(int i = 0; i < classNode.size(); i++) {
			classCount[i] = 0;
			for(int j = 0; j < features.get(0).size(); j++) {
				classCount[i] += trainResult[0][i][j];
			}
		}
		
		classDistrib.set(classCount);
		head.put(classDistrib);
		
		bn.add(head);
		List<CPT> featureList = new LinkedList<CPT>();
		for(int i = 0; i < windowSize; i++) {
			featureList.add(new CPT(features.get(i),classNode));
			for(int j = 0; j < classNode.size(); j++){
				EnumDistrib distrib = new EnumDistrib(Enumerable.aacid);
				distrib.set(trainResult[i][j]);
				featureList.get(i).put(distrib,classNode.getDomain().get(j));
			}
			bn.add(featureList.get(i));
		}
		
		VarElim ve = new VarElim();
        ve.instantiate(bn);
        Query q = null;
        double rightCase = 0;
        Iterator<Character[]> iterator = testing.iterator();
        Character[] record;
        while(iterator.hasNext()) {
			record = iterator.next();
			for(int i = 0; i < windowSize; i++) {
				featureList.get(i).setInstance(record[i]);
			}
			q = ve.makeQuery(classNode);
			CGTable r0 = (CGTable)ve.infer(q);
			double p = 0;
			int result = 0;
            for(int i = 0; i < classNode.size(); i++){
            	if(p < r0.getFactor(i)) {
            		result = i;
            		p = r0.getFactor(i);
            	}
            }
            if(classNode.getIndex(record[windowSize]) == result) {
            	rightCase ++;
            }
        }
		System.out.println(rightCase / testing.size());
	}
	
	
	
	public static double[][][] learn(List<Character[]> data, List<EnumVariable> features, EnumVariable parent) {
		double[][][] result = new double[features.size()][parent.size()][features.get(0).size()];
		for(int i = 0; i < features.size(); i++) {
			for(int j = 0; j < parent.size(); j++) {
				for(int x = 0; x < features.get(0).size(); x++) {
					result[i][j][x] = 0;
				}
			}
		}
		
		Iterator<Character[]> iterator = data.iterator();
		while(iterator.hasNext()) {
			Character[] record = iterator.next();
			int parentIdx = parent.getIndex(record[features.size()]);
			for(int i = 0; i < features.size(); i++) {
				int valueIdx = features.get(i).getIndex(record[i]);
				result[i][parentIdx][valueIdx] ++;
			}
		}
		
		return result;
	}
	
	/**
	 * load data from scv file
	 * @param fileName the data file
	 * @param attributeSize the attribute number
	 * @return the list of all attribute and label for class
	 */
	public static List<Character[]> loadData(String fileName, int attributeSize) {
		List<Character[]> resultList = new LinkedList<Character[]>();
		String line = "";
		BufferedReader br = null;
		try{
			br = new BufferedReader(new FileReader(fileName));
		
			while ((line = br.readLine()) != null) { 
				// use comma as separator
				String[] content = line.split(",");
				Character[] data  = new Character[attributeSize + 1];
				for(int i = 0; i < attributeSize; i++) {
					data[i] = content[0].charAt(i);
				}
				data[attributeSize] = content[1].charAt(0);
				resultList.add(data);
			}
			
		}catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (br != null) {
				try {
					br.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
		return resultList;
	}
	
}
