package bn.example;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Iterator;
import java.util.Map;

import dat.Continuous;
import dat.Domain;
import dat.DoubleSeq;
import dat.EnumVariable;
import dat.Enumerable;
import dat.IntegerSeq;
import dat.Variable;
import bn.prob.EnumDistrib;
import bn.alg.CGTable;
import bn.alg.EM;
import bn.alg.MAP;
import bn.alg.Query;
import bn.alg.VarElim;
import bn.node.CPT;
import bn.node.DirDT;
import bn.BNet;
import bn.Predef;
/**
 * used for secondary structure prediction
 * with algiment 
 * @author wangyufei
 */

public class SecondaryStructurePredictor {
	
	//final static private int windowSize = 13;
	final static private double trainDataRadio = 0.6f;
	//Naive Bayse
	final static private String fileName = "data/points.csv";
	static private List<String> testAmiodData;
	static private List<String> testSecondaryData;
	//Dirichlet Distribution
	final static private String dataSetFolder = "/Users/wangyufei/Desktop/100RESULT";
	final static private int attributesNum = 20;
	
	public static void main(String[] args) {
		/*
		for(int i = 5; i < 19; i+=2) {
			double result = 0.0;
			for(int j = 0; j < 3; j++) {
				result += DirPrediction(i);
			}
			System.out.println("window size = " + String.valueOf(i) + " " + String.valueOf(result / 3));
		}*/
		/*
		Enumerable enumerable = new Enumerable(new Character[]{'A', 'C', 'B', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'O', 'N', 'Q', 'P', 'S', 'R', 'U', 'T', 'W', 'V', 'Y', 'X', 'Z','-'});
		for(int i = 0; i < enumerable.size(); i++){
			System.out.println(enumerable.get(i));
		}*/
		
		
		SecondaryStructurePredictor p = new SecondaryStructurePredictor();
		p.start();
		
		
		//List<List<int[]>> data = loadMSAData();
		/* The data, the windowSize, and boolean value if output procession information*/
		//System.out.println(DirPrediction(data, 13, true));
	}
	
	public void start() {
		List<List<int[]>> data = loadMSAData();
		int size = data.size() / 10;
		sudocountForSequence(data);
		Collections.shuffle(data);
		DirDemo[] demos = new DirDemo[10];
		
		for(int i = 0; i < 10; i++) {
			List<List<int[]>> testing = data.subList(size * i, size * (i+1));
			List<List<int[]>> training = new LinkedList<List<int[]>>();
			for(List<int[]> a: data) {
				if(!testing.contains(a)) {
					training.add(a);
				}
			}
			String name = "runer" + String.valueOf(i);
			demos[i] = new DirDemo(name,training,testing, 13, false);
			demos[i].run();
		}
	}
	
	
	public static double NaiveBayes(int windowSize) {
		EnumVariable classNode = Predef.SecondaryStructure("secondary structure class");
		BNet bn = new BNet();
		
		//train
		List<EnumVariable> features = new LinkedList<EnumVariable>();
		initProcessData("data/result.out",windowSize);
		List<Character[]> training = loadData("data/points.csv",windowSize);
		
		for(int i = 0; i < windowSize; i++) {
			features.add(Predef.AminoAcid("feature" + String.valueOf(i)));
		}
		//double[][][] trainResult = learn(training,features,classNode);
		
		CPT head = new CPT(classNode);
		//EnumDistrib classDistrib = new EnumDistrib(classNode.getDomain());
		List<CPT> featureList = new LinkedList<CPT>();
		for(int i = 0; i < windowSize; i++) {
			featureList.add(new CPT(features.get(i),classNode));
		}
		
		
		Iterator<Character[]> iterator = training.iterator();
		//Character[] record;
		Character[][] datasetforEM = new Character[training.size()][];
		int count = 0;
		while(iterator.hasNext()) {
			//record = iterator.next();
			/*
			for(int i = 0; i < windowSize; i++) {
				featureList.get(i).countInstance(new Object[] {record[windowSize]}, record[i]);
			}
			head.countInstance(null, record[windowSize]);
			*/
			datasetforEM[count ++] = iterator.next();
		}
		
		long seed = 1;
		
		
		for(int i = 0; i < windowSize; i++) {
			//featureList.get(i).maximizeInstance();
			bn.add(featureList.get(i));
		}
		//head.maximizeInstance();
		bn.add(head);
		
		
		MAP em = new MAP(bn);
		Variable[] var = new Variable[windowSize + 1];
		for(int i = 0; i < windowSize; i++) {
			var[i] = features.get(i);
		}
		var[windowSize] = classNode;
        em.train(datasetforEM, var, seed);
		
		VarElim ve = new VarElim();
        ve.instantiate(bn);
        Query q = null;
        double rightCase = 0;
        double totalCase = 0;
        Iterator<String> iter1 = testAmiodData.iterator();
        Iterator<String> iter2 = testSecondaryData.iterator();
        String amiodData;
        String secondaryData;
        EnumDistrib distrib;
        while(iter1.hasNext() && iter2.hasNext()) {
        	amiodData = iter1.next();
        	secondaryData = iter2.next();
        	for(int i = 0; i < amiodData.length() - windowSize + 1; i++) {
        		for(int j = 0; j < featureList.size(); j++){
        			featureList.get(j).setInstance((Character)amiodData.charAt(i + j));
        		}
        		q = ve.makeQuery(classNode);
    			CGTable r0 = (CGTable)ve.infer(q);
    			distrib = (EnumDistrib)r0.query(classNode);
                totalCase ++;
                if(classNode.getIndex((Character)secondaryData.charAt(i + windowSize / 2)) 
                		== distrib.getMaxIndex()) {
                	rightCase ++;
                }
			}
			
        }
		File toDelete = new File(fileName);
		if(toDelete.exists() && toDelete.isFile()) {
			toDelete.delete();
		}
		return rightCase / totalCase;
	}
	
	private class DirDemo extends Thread {
		private int windowSize;
		private boolean isLog;
		private int attributesNum = 20;
		private List<List<int[]>> training;
		private List<List<int[]>> testing; 
		public DirDemo(String name, List<List<int[]>> trdata,List<List<int[]>> tedata, int windowSize, boolean isLog) {
			super(name);
			this.training = trdata;
			this.testing = tedata;
			this.windowSize = windowSize;
			this.isLog = isLog;
		}
		
		public void run() {
			System.out.println(getName() + " start !");
			// construct BN
			BNet bn = new BNet();
			
			EnumVariable classNode = Predef.SecondaryStructure("secondary structure class");
			List<Variable<EnumDistrib>> features = new LinkedList<Variable<EnumDistrib>>();
			for(int i = 0; i < windowSize; i++) {
				EnumDistrib enumDistrib = new EnumDistrib(new Enumerable(attributesNum));
				features.add(new Variable<EnumDistrib>(enumDistrib,"feature" + String.valueOf(i)));
			}
			CPT head = new CPT(classNode);
			List<DirDT> featureList = new LinkedList<DirDT>();
			for(int i = 0; i < windowSize; i++) {
				featureList.add(new DirDT(features.get(i),classNode));
			}
			
			int trainingSequence = training.size();
			int count = 0;
			int y = 0;
			int seCount = 0;
			//training
			Iterator<List<int[]>> iter = training.iterator();
			int trainingPoint = 0;
			for(List<int[]> protein: training) {
				trainingPoint += (protein.size() - windowSize + 1);
			}
			Object[][] training_EM = new Object[trainingPoint][windowSize + 1];
			Object[][] training_sec = new Object[trainingPoint][windowSize + 1];
			while(iter.hasNext()) {
				List<int[]> protein = iter.next();
				/*
				for(int[] record: protein) {
					Object secondary = classNode.getDomain().get((int)record[record.length - 1]);
					head.countInstance(null, secondary);
				}*/
				Object[] second = new Object[protein.size()];
				for(int i = 0; i < protein.size(); i++) {
					int[] classPoint = protein.get(i);
					second[i] = classNode.getDomain().get(classPoint[classPoint.length - 1]);
				}
				for(int i = 0; i <= protein.size() - windowSize; i++) {
					int[] classPoint = protein.get(i + windowSize / 2);
					training_EM[count][windowSize] 
							= classNode.getDomain().get(classPoint[classPoint.length - 1]);
					for(int j = 0; j < windowSize; j++) {
						int[] hist = protein.get(i + j);
						int[] attributes = new int[attributesNum];
						for(int x = 0; x < attributesNum; x++) {
							attributes[x] = hist[x];
						}
						training_EM[count][j] = IntegerSeq.intSeq(attributes);
					}
					count ++;
				}
				
				for(int i = 0; i <= second.length - windowSize; i++) {
					for(int j = 0; j < windowSize; j++) {
						training_sec[seCount][j] = second[i+j];
					}
					training_sec[seCount][windowSize] = second[i + windowSize / 2];
					seCount ++;
				}
				
				if(isLog){
					System.out.println("Finish training: " + String.valueOf(++y) + "/" + String.valueOf(trainingSequence));
				}
				
			}
			//head.maximizeInstance();
			bn.add(head);
			for(int i = 0; i < windowSize; i++) {
				//featureList.get(i).maximizeInstance();
				bn.add(featureList.get(i));
				/*
				if(isLog){
					System.out.println("Finish node training: " + String.valueOf(i + 1) + "/" + String.valueOf(windowSize));
				}*/
			}
			
			//layer two
			BNet bn2 = new BNet();
			MAP em2 = new MAP(bn2);
			EnumVariable layerTwoClassNode = Predef.SecondaryStructure("secondary structure class");
			CPT predictor = new CPT(layerTwoClassNode);
			List<EnumVariable> sfeatures = new LinkedList<EnumVariable>();
			List<CPT> sfeatureList = new LinkedList<CPT>();
			Variable[] vars2 = new Variable[windowSize + 1];
			for(int i = 0; i < windowSize; i++) {
				sfeatures.add(Predef.SecondaryStructure("feature" + String.valueOf(i)));
				sfeatureList.add(new CPT(sfeatures.get(i),layerTwoClassNode));
				bn2.add(sfeatureList.get(i));
				vars2[i] = sfeatures.get(i);
			}
			vars2[windowSize] = layerTwoClassNode;
			bn2.add(predictor);
			em2.train(training_sec, vars2, 1);
			VarElim ve2 = new VarElim();
	        ve2.instantiate(bn2);
			
			MAP em = new MAP(bn);
			Variable[] var = new Variable[windowSize + 1];
			for(int i = 0; i < windowSize; i++) {
				var[i] = features.get(i);
			}
			var[windowSize] = classNode;
	        em.train(training_EM, var, 1);
			
			if(isLog){
				System.out.println("Start to test");
			}
			int testingSequence = testing.size();
			int tested = 0;
			//testing
			VarElim ve = new VarElim();
	        ve.instantiate(bn);
	        Query q = null;
			double totalCase = 0;
			double rightCase = 0;
			EnumDistrib distrib;
			iter = testing.iterator();
			while(iter.hasNext()) {
				List<int[]> protein = iter.next();
				Object[] result = new Object[protein.size()];
				for(int i = 0; i <= protein.size() - windowSize; i++) {
					for(int j = 0; j < windowSize; j++) {
						int[] hist = protein.get(i + j);
						int[] attributes = new int[attributesNum];
						for(int x = 0; x < attributesNum; x++) {
							attributes[x] = hist[x];
						}
						featureList.get(j).setInstance(IntegerSeq.intSeq(attributes));
					}
					q = ve.makeQuery(classNode);
	    			CGTable r0 = (CGTable)ve.infer(q);
	    			distrib = (EnumDistrib)r0.query(classNode);
	    			result[i + windowSize / 2] = distrib.getMax();
				}
				for(int x = 0; x < windowSize / 2; x++) {
					result[x] = result[windowSize / 2];
					result[protein.size() - 1 - x] = result[protein.size() - windowSize / 2 - 1];
				}
				for(int i = 0; i <= result.length - windowSize; i++) {
					for(int j = 0; j < windowSize; j++) {
						sfeatureList.get(j).setInstance(result[i + j]);
					}
					q = ve2.makeQuery(layerTwoClassNode);
					CGTable r0 = (CGTable)ve2.infer(q);
	    			distrib = (EnumDistrib)r0.query(layerTwoClassNode);
	    			if(distrib.getMax().equals(result[i + windowSize / 2])) {
	    				rightCase ++;
	    			}
	    			totalCase ++;
				}
				if(isLog) {
					System.out.println("Finish training: " + String.valueOf(++tested) + "/" + String.valueOf(testingSequence));
				}
			}
			
			System.out.println(getName() + "  result : " + String.valueOf(rightCase / totalCase));
		}
	}
	
	
	
	public static void sudocountForSequence(List<List<int[]>> data) {
		Iterator<List<int[]>> iter = data.iterator();
		Iterator<int[]> iter1;
		while(iter.hasNext()) {
			List<int[]> protein = iter.next();
			iter1 = protein.iterator();
			while(iter1.hasNext()) {
				int[] record = iter1.next();
				for(int i = 0; i < record.length - 1; i++) {
					//record[i] *= 100;
					record[i] += 1;
				}
			}
		}
	}
	
	public static List<List<int[]>> loadMSAData() {
		List<List<int[]>> result = new LinkedList<List<int[]>>();
		File folder = new File(dataSetFolder);
		for(File file: folder.listFiles()) {
			List<int[]> list = loadSingleMSAData(file.getAbsolutePath());
			result.add(list);
		}
		return result;
	} 
	
	public static List<int[]> loadSingleMSAData(String fileName)  {
		List<int[]> result = new LinkedList<int[]>();
		BufferedReader br = null;
		String line;
		String[] stringData;
		EnumVariable classNode = Predef.SecondaryStructure("secondary structure class");
		try {
			br = new BufferedReader(new FileReader(fileName));
			while ((line = br.readLine()) != null) {
				stringData = line.split(",");
				int[] data = new int[attributesNum + 1];
				for(int i = 0; i < attributesNum; i++) {
					data[i] = Integer.valueOf(stringData[i]);
				}
				
				Character s = stringData[stringData.length - 1].charAt(0);
				data[attributesNum] = classNode.getDomain().getIndex(s);
				result.add(data);
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return result;
	}
	
	public static void initProcessData(String initDataFileName, int windowSize) {
		BufferedReader br = null;
		BufferedWriter bw = null;
		String amiod = "";
		String secondary = "";
		String[] data;
		List<String> initData = new LinkedList<String>();
		List<String> initTrainingData ;
		List<String> initTestData;
		try{
			br = new BufferedReader(new FileReader(initDataFileName));
			bw = new BufferedWriter(new FileWriter(new File(fileName))); 
			while ((amiod = br.readLine()) != null) { 
				secondary = br.readLine();
				initData.add(amiod + ";" + secondary);
			}
			Collections.shuffle(initData);
			int cutOffIndex = (int)(initData.size() * trainDataRadio);
			initTrainingData = initData.subList(0, cutOffIndex);
			initTestData = initData.subList(cutOffIndex, initData.size());
			
			Iterator<String> iter = initTrainingData.iterator();
			 
			while(iter.hasNext()) {
				data = iter.next().split(";");
				amiod = data[0];
				secondary = data[1];
				
				for(int i = 0; i < amiod.length() - windowSize + 1; i++) {
					bw.write(amiod.substring(i, i + windowSize) + "," + secondary.charAt(i + windowSize / 2));
					bw.newLine();  
	                bw.flush();  
				}
			}
			
			testAmiodData = new LinkedList<String>();
			testSecondaryData = new LinkedList<String>();
			
			iter = initTestData.iterator();
			while(iter.hasNext()) {
				data = iter.next().split(";");
				testAmiodData.add(data[0]);
				testSecondaryData.add(data[1]);
			}
			
		}catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (br != null) {
				try {
					br.close();
					bw.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	/*
	 * not necessary
	public static double[][][] learn(List<Character[]> data, List<EnumVariable> features, EnumVariable parent) {
		double[][][] result = new double[features.size()][parent.size()][];
		for(int i = 0; i < features.size(); i++) {
			for(int j = 0; j < parent.size(); j++) {
				result[i][j] = new double[features.get(i).size()];
				for(int x = 0; x < features.get(i).size(); x++) {
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
	*/
	
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
