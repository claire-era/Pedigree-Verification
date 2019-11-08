package pedigreeVerification;

import java.util.ArrayList;

import org.rosuda.JRI.Rengine;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
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
		int hetCount, site = 0;
		int numOfSites = genos.numberOfSites(); // for 0-indexing
		int[] hcArr = new int[numOfSites];

		for (site = 0; site < numOfSites; site++) {
			hetCount = genos.heterozygousCount(site);
			hcArr[site] = hetCount;
		}

		// get indices of max count of polymorphism
		int max = hcArr[0];
		for (int i = 1; i < hcArr.length; i++)
			if (hcArr[i] > max)
				max = hcArr[i];
		

		// find all indices of polymorphic markers matching the max value
		for (site = 0; site < numOfSites; site++)
			if (genos.isPolymorphic(site) && max == genos.heterozygousCount(site)) {
				polyMarkers.add(site);
			}
		

		return polyMarkers; // returns index of the location of the polymorphic marker
	}

	private static ArrayList<Integer> getInitSolution(ArrayList<Integer> polymarkers) {
		ArrayList<Integer> taxa = new ArrayList<Integer>();

		int psite = polymarkers.get(0); // get the max site
//		System.out.println(psite + " psite to ");
		for (int taxon = 0; taxon < genos.numberOfTaxa(); taxon++)
			if (genos.isHeterozygous(taxon, psite)) {
//				System.out.println(taxon);
				taxa.add(taxon); // adds index of a sample included in the initial solution
			}

		return taxa;
	}

	private static String createVector(double[] observed, double[] expected) {
		// CREATE VECTOR OF (OBSERVED,EXPECTED) values
		String res = "c(";
		for (int a = 0; a < observed.length; a++)
			res = res + Double.toString(observed[a]) + ",";

		for (int a = 0; a < expected.length; a++) {
			if (a == (expected.length - 1))
				res = res + Double.toString(expected[a]) + ")";
			else
				res = res + Double.toString(expected[a]) + ",";
		}
		// -----------------
		return res;
	}

	private static Double getPvalue(Rengine re, double[] observed, double[] expected) {
		String res = createVector(observed, expected);
		Double pval = null;
		if (!re.waitForR()) {
			System.out.println("Cannot load R");
			return null;
		}
		try {
			// ----- evaluate string
			String sample_m = "matrix(" + res + ",nrow = 3)"; //for each element inside the classification array. expected[0]..[2].
			re.eval("observed=" + sample_m);
			re.eval("result=fisher.test(observed, alternative = \"greater\")");
			pval = re.eval("result$p.value").asDouble();
		} catch (Exception e) {
		}
		return pval;
	}
	
	private static int getMax(double[] arr) { //returns index of max value in an array
		double max = arr[0];
		int i = 0;
		int max_i = 0;
		for (i = 1; i < arr.length-1; i++)
			if (arr[i] > max) {
				max = arr[i];
				max_i = i;
			}
//		System.out.println(i);
		return max_i;

	}
	
	private static double countAlleles(ArrayList<Integer> listOftaxa, int site, GenotypeTable genos, int type) {
		double majCount = 0;
		double minCount = 0;
		double hetCount = 0;
		byte majorAlleleInSite = genos.majorAllele(site);
		byte minorAlleleInSite = genos.minorAllele(site);
		
		for(int i = 0; i < listOftaxa.size(); i++) {
			int t = listOftaxa.get(i);
			if(genos.isHeterozygous(t, site)) hetCount++;
			else if(majorAlleleInSite == genos.genotypeArray(t,site)[1])	majCount++;
			else if(minorAlleleInSite == genos.genotypeArray(t,site)[1]) minCount++;
		}
		if(type == 0) return majCount;
		else if (type == 1) return minCount;
		else return hetCount;
	}

	private static ArrayList<Integer> GetCrossType(GenotypeTable genos) {
		ArrayList<Integer> ct = new ArrayList<Integer>(); // CROSS TYPE / ALLELE COMBINATION after function
		ArrayList<Integer> polymarkers = GetPolymorphicMarkers(genos);
		ArrayList<Integer> initTaxaSol = getInitSolution(polymarkers); // generated initial solution after good marker
																		// selection //taxa indices of initial solution
		Rengine re = new Rengine(null, false, null);
		int nCountSites = genos.numberOfSites();
		double initSolTaxa = initTaxaSol.size();
		double X, Y, Z;
		
		double[] expected1 = {initSolTaxa, 0, 0};
		double[] expected2 = {0, 0, initSolTaxa};
		double[] expected3 = {0, initSolTaxa, 0 }; // ASSUMPTION: ALL SAMPLES ARE HET (Segregation: AA x aa) 100% het
		double[] expected4 = {initSolTaxa / 2, initSolTaxa / 2, 0 }; // ASSUMPTION: HALF ARE HET AND HALF ARE HOMO	// DOMINANT (Segregation: AA x Aa) 50% self
		double[] expected5 = {0, initSolTaxa/2, initSolTaxa/2};
		double[] expected6 = {initSolTaxa/4, initSolTaxa/2, initSolTaxa/4};
		
		
		for (int site = 0; site < nCountSites; site++) {
			X = countAlleles(initTaxaSol, site, genos, 0); //0 for major allele
			Z = countAlleles(initTaxaSol, site, genos, 1); //1 for minor allele
			Y = countAlleles(initTaxaSol, site, genos, 2); // HETEROZYGOUS COUNT
			
			double[] observed = { X, Y, Z };
			
			double p_t1 = Step1.getPvalue(re, observed, expected1);
			double p_t2 = Step1.getPvalue(re, observed, expected2);
			double p_t3 = Step1.getPvalue(re, observed, expected3);
			double p_t4 = Step1.getPvalue(re, observed, expected4);
			double p_t5 = Step1.getPvalue(re, observed, expected5);
			double p_t6 = Step1.getPvalue(re, observed, expected6);
			
			double[] p_arr = {p_t1, p_t2, p_t3, p_t4, p_t5, p_t6};
			
			//classify: get max between 6 numbers
			int index_max = getMax(p_arr);
			index_max = index_max + 1; //for fixing 0-indexing, classification
			
			ct.add(index_max);
		}
		
		re.end(); //END R ENGINE INSTANCE AFTER PERFORMING fisher.test() ON ALL SITES/SNPS
		return ct;
	}
	
	private static void InferGenotype(GenotypeTable genos, ArrayList<Integer> crossType) {
		// for every case, there must be some metric/standard/expectation to be followed
		
		//but FIRST: set STANDARD
		byte[] exp_f1 = new byte[2]; //exp0 holds higher/major bit,1 holds lower/min/het
		byte[] exp_self = new byte[2];
		byte majAllele;
		byte minAllele;
		for(int site = 0; site < crossType.size(); site++) { //for each SNP, with corresponding crosstype
			int cross = crossType.get(site); //get crosstype per site
			majAllele = genos.majorAllele(site);
			minAllele = genos.minorAllele(site);
			//INCUR DOMINANCE / majority / minority
			if(cross == 1) {  //MAJ X MAJ
				exp_f1[0] = majAllele;
				exp_f1[1] = majAllele;
				
				exp_self[0] = ;
			}else if(cross == 2) { //MIN X MIN
				exp_f1[0] = minAllele;
				exp_f1[1] = minAllele;
			}else if(cross == 3) { // MAJ X MIN == HET
				exp_f1[0] = majAllele;
				exp_f1[1] = minAllele;
			}else if(cross == 4) { //MAJOR x HET
				exp_f1[0] = majAllele;
				exp_f1[1] = GenotypeTableUtils.getUnphasedDiploidValue(majAllele, minAllele);
//				System.out.println(NucleotideAlignmentConstants.getNucleotideIUPAC(exp_f1[1]));
 			}else if(cross == 5) { //HET X MINOR
				exp_f1[0] = GenotypeTableUtils.getUnphasedDiploidValue(majAllele, minAllele);
				exp_f1[1] = minAllele;
			}else if(cross == 6) { //HET X HET
				exp_f1[0] = GenotypeTableUtils.getUnphasedDiploidValue(majAllele, minAllele);
				exp_f1[1] = exp_f1[0];
			}
			
//			System.out.println(genos.genotypeArray(1, 1));
		}
		
	}

	private static void checkInputs(String[] args) {
		try {
			if (args.length == 0) {
				System.out.println(
						"No arguments. \nPlease supply .hmp.txt file and pedigree file separated by space.\nFormat: java -jar PedVer.jar <input.hmp.txt> <pedigree_file.txt> <cut-off>");
			} else if (args.length == 3) { // then there must be a provided cut-off
				start(args[0], args[1], Double.parseDouble(args[2]));
			} else if (args.length == 2) { // then default, get minimum
				start(args[0], args[1], 0.0);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void ProgramEnded() {
//		System.out.println("... Done.");
//		Thread.sleep(500);
		long startTime = System.nanoTime();
		long endTime = System.nanoTime();
		System.out.println("Took " + (endTime - startTime) + " ns");
	}

	private static void start(String hmpFile, String pedFile, double cutOff) throws InterruptedException {
		System.out.println("Cut-off: " + cutOff);

		try {
			// READ FILE
//			Step1.listOfGenos = ReadFile(hmpFile, pedFile);
			// GET KEYS FROM CORRESPONDING PEDIGREE FILE
//			ArrayList<String> keys = PedigreeFileInfo.getKeys();

			// IMPLEMENT STEP 1
			genos = Step1.ReadHMPFile(hmpFile);
			ArrayList<Integer> crossType = Step1.GetCrossType(genos);
			Step1.InferGenotype(genos, crossType);
			
			for(int i = 0; i < crossType.size(); i++) System.out.println(crossType.get(i));

			ProgramEnded();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
