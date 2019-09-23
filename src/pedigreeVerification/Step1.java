package pedigreeVerification;

import java.io.IOException;
import java.util.ArrayList;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;


public class Step1 {
	private static ArrayList<GenotypeTable> listOfGenos;
	private static GenotypeTable genos;
	
	public Step1(String[] args) {
		Step1.checkInputs(args);
	}
	
	private static ArrayList<GenotypeTable> ReadFile(String hmpFile, String pedFile) { // Read file contents
		ArrayList<GenotypeTable> listOfGenos = PedigreeFileInfo.getFilteredGenotypeTable(hmpFile, pedFile);
		return listOfGenos;
	}
	
	private static GenotypeTable ReadHMPFile(String hmp_filename) {
		genos = ImportUtils.readFromHapmap(hmp_filename);
		return genos;
	}

	private static ArrayList<Integer> GetPolymorphicMarkers(GenotypeTable genos) {
		ArrayList<Integer> polyMarkers = new ArrayList<>();

		int numOfSites = genos.numberOfSites(); // for 0-indexing

		for (int site = 0; site < numOfSites; site++) {
			if (genos.isPolymorphic(site))	polyMarkers.add(site);
		}
		return polyMarkers; //returns index of the location of the polymorphic marker
	}
	
	private static ArrayList<Integer> GetCrossType(GenotypeTable genos){
		ArrayList<Integer> ct = new ArrayList<Integer>();
		ArrayList<Integer> polymarkers = GetPolymorphicMarkers(genos);
		int psite;
		double p, q, X, Y, Z;
		int N;
		double[] observed = new double[3];
		double[] expected1 = {0,1,0};
		double[] expected2 = {0.5,0.5,0};
		double[][] result;
		
		for(int i = 0; i < polymarkers.size(); i++) {
			psite = polymarkers.get(i);
			p = genos.majorAlleleFrequency(psite);
			q = genos.minorAlleleFrequency(psite);
			N = genos.numberOfTaxa();
			Y = genos.heterozygousCount(psite);
			X = ((2*N*p) - Y)/2;
			Z = ((2*N*q) - Y)/2;
			X = X/N;
			Z = Z/N;
			Y = Y/N;
			
			//CLASSIFY IF: TYPE 1 or 2
			//APPLY FISHER'S TEST to see if observed values correspond with the exact values
			observed[0] = X;
			observed[1] = Y;
			observed[2] = Z;
			
//			System.out.println(observed[0]+" : "+observed[1]+" : "+observed[2]);
//			result = ChiSquare(expected1, observed);
//			FisherExact fishers_test = new FisherExact(N);
//			System.out.println(fishers_test.getP(10, 0, 5, 4));
 		}
		

		return null;
	}
	
	private static void checkInputs(String[] args) {
		try {
			if(args.length==0) {
				System.out.println("No arguments. \nPlease supply .hmp.txt file and pedigree file separated by space.\nFormat: java -jar PedVer.jar <input.hmp.txt> <pedigree_file.txt> <cut-off>");
			}else if(args.length==3){ //then there must be a provided cut-off
				start(args[0], args[1], Double.parseDouble(args[2]));
			}else if(args.length==2) { //then default, get minimum
				start(args[0],args[1],0.0);
			}
		}catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void ProgramEnded() {
//		System.out.println("... Done.");
//		Thread.sleep(500);
		long startTime = System.nanoTime();
		long endTime = System.nanoTime();
		System.out.println("Took "+(endTime - startTime) + " ns");
	}
	
	private static void start(String hmpFile, String pedFile, double cutOff) throws InterruptedException {
		System.out.println("Cut-off: " +  cutOff);
		
		try {
			//READ FILE
//			Step1.listOfGenos = ReadFile(hmpFile, pedFile);
			//GET KEYS FROM CORRESPONDING PEDIGREE FILE
//			ArrayList<String> keys = PedigreeFileInfo.getKeys();
			
			//IMPLEMENT STEP 1
			genos = Step1.ReadHMPFile(hmpFile);
			Step1.GetCrossType(genos);
			
			ProgramEnded();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
	

