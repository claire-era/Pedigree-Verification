package pedigreeVerification;

import java.io.IOException;
import java.util.ArrayList;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;

public class PedVerification {
	public static ArrayList<GenotypeTable> ReadFile() { // Read file contents
		ArrayList<GenotypeTable> listOfGenos = PedigreeFileInfo.getFilteredGenotypeTable();
		return listOfGenos;
	}

	private static ArrayList<Integer> GetPolymorphicMarkers(GenotypeTable genos) {
		ArrayList<Integer> polyMarkers = new ArrayList<>();

		int numOfSites = genos.numberOfSites(); // for 0-indexing

		for (int site = 0; site < numOfSites; site++) {
			if (genos.isPolymorphic(site)) {
				polyMarkers.add(site);
			}
		}
		

		return polyMarkers;
	}

	private static byte[] GetCrossType(GenotypeTable genos, ArrayList<Integer> polyMarkers) throws IOException {
		byte[] snp;
		int polyMarkersLength = polyMarkers.size();
		byte[] crossTypesAllSites = new byte[polyMarkersLength];
		int site;
		int numOfTaxa = genos.numberOfTaxa();
		long hetC = 0, homoDC = 0, homoRC = 0, nanC = 0; // counters for genotypes
		byte hetG = 0, homoDG = 0, homoRG = 0;
		String majAllele, minAllele;
		String geno;
		double[] observedVal = new double[3];

		double[] expectedValType1, expectedValType2, expectedValType3;

		for (int taxa = 0; taxa < polyMarkersLength; taxa++) { // for every SNP, check genotypic frequency
			hetC = homoDC = homoRC = 0;
			site = polyMarkers.get(taxa);
			snp = genos.genotypeAllTaxa(site);
//			System.out.println(genos.siteName(site));

			// Obtain all heterozygous genotypes
			hetC = genos.heterozygousCount(site);

			majAllele = genos.majorAlleleAsString(site);
			minAllele = genos.minorAlleleAsString(site);

			// Obtain all homozygous genotypes
			for (int u = 0; u < snp.length; u++) {
				geno = genos.genotypeAsString(u, site);
				if (snp[u] == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
					nanC++;
					continue;
				}
				if (geno.toLowerCase().contains(majAllele.toLowerCase())) { // if gene is dominant
					homoDG = genos.genotype(u, site);
					homoDC++;
				}
				if (geno.toLowerCase().contains(minAllele.toLowerCase())) { // if gene is recessive
					homoRG = genos.genotype(u, site);
					homoRC++;
				}
				if (genos.isHeterozygous(u, site))
					hetG = genos.genotype(u, site);
			}

			expectedValType1 = calculateExpectedValues(numOfTaxa, 1);
			expectedValType2 = calculateExpectedValues(numOfTaxa, 2);
			expectedValType3 = calculateExpectedValues(numOfTaxa, 3);
			observedVal[0] = homoDC;
			observedVal[1] = hetC;
			observedVal[2] = homoRC;

			// VERIFY PROBABILITIES AND IDENTIFY GENOTYPE

			double[][] expected = new double[3][3];

			expected[0] = expectedValType1;
			expected[1] = expectedValType2;
			expected[2] = expectedValType3;

			crossTypesAllSites[taxa] = getCross(expected, observedVal);
		}

		return crossTypesAllSites;
	}

	private static byte getCross(double[][] expected, double[] observed) {

		double[] result = FishersTest(expected, observed); // returns [.24,.67,.09] for [type1, type2, type3]
		double highestP = result[0];
		byte cross = 0;
		if (result[1] > highestP) {
			highestP = result[1];
			cross = 1;

		}
		if (result[2] > highestP) {
			highestP = result[2];
			cross = 2;
		}
		return cross;
	}

	private static double[] FishersTest(double[][] expected, double[] observed) {
		double[] result = null;// = new double[expected.length]; //3 cross types
//		Random r = new Random();
		// 0: 1: 2::Type 1:Type 2: Type 3

//		for(int i = 0; i < result.length; i++){
//			result[i] = ChiSquare(expected, observed, i); // i - type  ////////////////////// fix!! iterate thru expected[][] to get actual values
//			result[i] = r.nextDouble();

//			System.out.println(result[i]);
//		}

		for (int i = 0; i < expected.length; i++) {
//			result[i] = ChiSquare(expected[i], observed, i);
//			if(i == 0 && result==null) { //type 1
//		result = new double[] { 0.90, 0.50, 0.30 };
//				continue;
//			}if(i == 1 && result==null) { //type 2
//		result = new double[] { 0.50, 0.90, 0.30 };
//				continue;
//			}if(i == 2 && result==null) { //type 3
//		result = new double[] { 0.30, 0.50, 0.90 };
//				continue;
//			}
//			System.out.println(observed[1] + " " + observed[0]);
			if (observed[1] >= observed[0]) {
				// type 1
				result = new double[] { 0.90, 0.50, 0.30 };
			} else if (observed[0] >= observed[1]) {
				// type 2
				result = new double[] { 0.50, 0.90, 0.30 };
			} else {
				// type 3
				result = new double[] { 0.30, 0.50, 0.90 };
			}
//			System.out.println(i);
		}

//		for (int i = 0; i < result.length; i++) {
//			result[i] = ChiSquare(expected, observed, i); // i - type  ////////////////////// fix!! iterate thru expected[][] to get actual values
//			result[i] = r.nextDouble();

//			System.out.println(result[i]);
//		}

		return result;
	}

	private static byte[][] getParentGenotypeAllSites(GenotypeTable genos, ArrayList<Integer> polyMarkers,
			byte[] crossAllSites) {
		byte[][] parentGenotypeAllSites = new byte[polyMarkers.size()][2]; // even length of array because first two
																			// elements correspond to a gene
		String p1_geno = null, p2_geno = null;
		String majA, minA;

		for (int i = 0; i < polyMarkers.size(); i++) {
			majA = genos.majorAlleleAsString(polyMarkers.get(i));
			minA = genos.minorAlleleAsString(polyMarkers.get(i));
			if (crossAllSites[i] == 0) { // type 1
				// AA x aa
				p1_geno = majA.concat(majA); // AA
				p2_geno = minA.concat(minA); // aa
			}

			else if (crossAllSites[i] == 1) { // type 2
				// AA x Aa
				p1_geno = majA.concat(majA); // AA
				p2_geno = majA.concat(minA); // Aa
			}

			else if (crossAllSites[i] == 2) { // type 3
				// Aa x Aa
				p1_geno = majA.concat(minA); // Aa
				p2_geno = majA.concat(minA); // Aa
			}

//			System.out.print(NucleotideAlignmentConstants.NUCLEOTIDE_IUPAC_HASH.get(NucleotideAlignmentConstants.getNucleotideDiploidByte(p1_geno)));
//			System.out.println(NucleotideAlignmentConstants.NUCLEOTIDE_IUPAC_HASH.get(NucleotideAlignmentConstants.getNucleotideDiploidByte(p2_geno)));
			parentGenotypeAllSites[i][0] = NucleotideAlignmentConstants.getNucleotideDiploidByte(p1_geno);
			parentGenotypeAllSites[i][1] = NucleotideAlignmentConstants.getNucleotideDiploidByte(p2_geno);
		}

		return parentGenotypeAllSites;
	}

	private static double[] calculateExpectedValues(int numOfTaxa, int typeNo) {
		double[] expectedValuesOfType = new double[3];
		// GIVEN THAT A:B:C corresponds to: AA:Aa:aa
		for (int i = 0; i < 3; i++) {
			if (typeNo == 1) { // follow 0:1:0 ratio
				expectedValuesOfType[0] = 0;
				expectedValuesOfType[1] = numOfTaxa;
				expectedValuesOfType[2] = 0;
			} else if (typeNo == 2) { // follow 0:0.5:0.5 ratio
				expectedValuesOfType[0] = numOfTaxa * 0.5;
				expectedValuesOfType[1] = numOfTaxa * 0.5;
				expectedValuesOfType[2] = 0;
			} else if (typeNo == 3) { // follow 0.25:0.5:0.25 ratio
				expectedValuesOfType[0] = numOfTaxa * 0.25;
				expectedValuesOfType[1] = numOfTaxa * 0.5;
				expectedValuesOfType[2] = numOfTaxa * 0.25;
			}

		}
		return expectedValuesOfType;
	}

	private static double[][] getScoresAllTaxa(GenotypeTable genos, ArrayList<Integer> polyMarkers,
			byte[] crossAllSites) {
		int polyMarkersLen = polyMarkers.size();
		double[][] scores = new double[polyMarkersLen][];
		int site = 0;
		byte[] snp;
		byte geno;
		String genoo;
		for (int s = 0; s < polyMarkersLen; s++) { // for every snp
			snp = genos.genotypeAllTaxa(s);
			scores[s] = new double[snp.length];
			site = polyMarkers.get(s);
//			System.out.println(snp.length);
			for (int taxon = 0; taxon < snp.length; taxon++) { // for every taxa/sample in the snp
//				System.out.println(genos.genotypeAsString(0, 1) + " " + NucleotideAlignmentConstants.NUCLEOTIDE_IUPAC_HASH.get(parentGenotypeAllSites[0][0]));		
//				System.out.println(genos.genotype(0, 1) + " " + parentGenotypeAllSites[0][0]);
				geno = genos.genotype(taxon, site);
				if (crossAllSites[s] == 0) { // if type 1
					if (genos.isHeterozygous(taxon, site)) {
						scores[s][taxon] = 0f;
					} else {
						scores[s][taxon] = 1.0f;
					}
				} else if (crossAllSites[s] == 1) { // if type 2
					if (genos.isHeterozygous(taxon, site)) { // het
						scores[s][taxon] = 0.05f;
					} else if (NucleotideAlignmentConstants.isHomozygousACGT(geno)) {
						genoo = genos.genotypeAsString(taxon, site);
						if (genoo.toLowerCase().contains(genos.majorAlleleAsString(s).toLowerCase())) { // dominant gene
							scores[s][taxon] = 0.1f;
						} else if (genoo.toLowerCase().contains(genos.minorAlleleAsString(s).toLowerCase())) { // recessive
																												// gene
							scores[s][taxon] = 1.0f;
						}
					}
				}
				/*
				 * else if (crossAllSites[s] == 2) { // if type 3 continue; }
				 */
			}

		}
		return scores;
	}

	private static double[] getMeanScoresAllTaxa(GenotypeTable genos, ArrayList<Integer> polyMarkers, double[][] scores) {
		double sum = 0;
		int i, j = 0;
		int rows = scores.length;
		int cols = scores[0].length;

		double[] scoresAllTaxa = new double[genos.numberOfTaxa()];
		for (j = 0; j < cols; j++) {
			for (i = 0; i < rows; i++) {
				sum += scores[i][j];
			}
			sum = sum / rows;
			scoresAllTaxa[j] = sum;
			sum = 0;
		}

		return scoresAllTaxa;
	}

	private static ArrayList<Integer> getOutCrossedTaxa(double[] meanScores) {// goal: get index of the minimum taxa
		ArrayList<Integer> indices = new ArrayList<Integer>();
		int minIndex = 0;

		// get minimum scores
		for (int i = 0; i < meanScores.length; i++) {
			if (meanScores[i] < meanScores[minIndex]) {
//				if (minIndex == 0 && ) {
//					indices.remove(minIndex);
//				}
				minIndex = i;
				indices.add(minIndex);
			} else if (meanScores[i] == meanScores[minIndex]) {
				indices.add(i);
			}
		}

		return indices;
	}

	private static ArrayList<Integer> getOutCrossedTaxa(double[] meanScores, double cutOff) {// goal: get index of the
																							// minimum taxa given that
																							// there is an indicated
																							// cutOff
		double lim = (double) cutOff;
		ArrayList<Integer> indices = new ArrayList<Integer>();
		byte minIndex = 0;

		// get minimum scores
		for (int i = 0; i < meanScores.length; i++) {
			System.out.println(meanScores[i]);
			if (meanScores[i] < meanScores[minIndex]) {
				if (meanScores[i] < lim)
					continue;
				if (indices.contains((int) minIndex)) {
					indices.remove(minIndex);
				}
				minIndex = (byte) i;
				indices.add(i);
			} else if (meanScores[i] == meanScores[minIndex]) {
				indices.add(i);
			}

		}

		return indices;
	}

	private static boolean[] verifyParentProgenyGenotype(GenotypeTable genos, ArrayList<Integer> outCrossTaxaIndex) {
		
		boolean[] accepted = new boolean[outCrossTaxaIndex.size()];
		/* for every accepted F1
		 * 		p1 =
		 * */
		PedigreeFileInfo.progenyGroupToParentMapping(genos, outCrossTaxaIndex);
		
		return accepted;
	}

	private static void print(GenotypeTable genos, ArrayList<Integer> polyMarkers, double[][] scores, double[] meanScores,
			ArrayList<Integer> outCrossTaxaIndex) {
//		System.out.print("\t");
//		for(int d = 0; d < genos.numberOfTaxa(); d++) {
//			System.out.print(genos.taxaName(d) + " ");
//		}
//		System.out.print("\n");
//		
//
//		for (int a = 0; a < scores.length; a++) {
//			System.out.print(genos.siteName(polyMarkers.get(a)) + ": ");
//			for (int b = 0; b < scores[a].length; b++) {
//				System.out.print(scores[a][b] + " ");
//			}
//			System.out.print("\n");
//		}
//		
//		System.out.print("\n");
//
//		for (int c = 0; c < meanScores.length; c++) {
//			System.out.println(genos.taxaName(c) + ": " + meanScores[c]);
//		}
//
		System.out.print("\n");

		System.out.println("Out-Crossed: ");
		for(int i = 0; i < outCrossTaxaIndex.size(); i++) {
			System.out.println(genos.taxaName(outCrossTaxaIndex.get(i)));
		}

	}

	public static void main(String[] args) {
		// STEP 2 Get polymorphic markers
		// STEP 3 Identify cross type: 1 or 2 or 3
		// For every marker inside polyMarkers, identify if Type 1 or 2 or 3 cross
		// Identify parent genotypes per site
		// Match taxa/sample according to the type of cross of the site and give score

		try {
			ArrayList<GenotypeTable> listOfGenos = ReadFile();
//			GenotypeTable parentGenos = PedigreeFileInfo.getParentsGenoTable();
			for (int i = 0; i < listOfGenos.size(); i++) {
				GenotypeTable genos = listOfGenos.get(i);
				ArrayList<Integer> polyMarkers = GetPolymorphicMarkers(genos);
				if (polyMarkers.size() != 0) {
					byte[] crossTypeAllSites = GetCrossType(genos, polyMarkers);
					byte[][] parentGenotypeAllSites = getParentGenotypeAllSites(genos, polyMarkers, crossTypeAllSites);
					double[][] scores = getScoresAllTaxa(genos, polyMarkers, crossTypeAllSites);
					double[] meanScores = getMeanScoresAllTaxa(genos, polyMarkers, scores);
					ArrayList<Integer> outCrossTaxaIndex = getOutCrossedTaxa(meanScores); // returns array of taxa that
																							// is the result of
																							// out-crossed
						
//					print(genos, polyMarkers, scores, meanScores, outCrossTaxaIndex);
					
					// step 2
					// boolean 2d array of every accepted F1, (or taxa indices of outCrossTaxaIndex)
					boolean[] acceptedTaxa = verifyParentProgenyGenotype(genos, outCrossTaxaIndex);
					break;
					
				}	
			}
			

			long startTime = System.nanoTime();
			long endTime = System.nanoTime();
			System.out.println("Took "+(endTime - startTime) + " ns"); 
			
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
