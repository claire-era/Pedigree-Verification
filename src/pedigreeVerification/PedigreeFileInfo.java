package pedigreeVerification;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.taxa.TaxaList;
/* A text file consisting of genotypic information.
 * with .hmp.txt **
 * */

public class PedigreeFileInfo {
	private static GenotypeTable genos;
	private static ArrayList<ArrayList<String>> taxaInfo;
	private static HashMap<String, ArrayList<Integer>> taxaPerGroup;

	public static int[][] progenyGroupToParentMapping(GenotypeTable currGenos, ArrayList<Integer> outCrossTaxaIndex) {
		String p1name;
		String p2name;
		int[][] progenyParentMap = new int[outCrossTaxaIndex.size()][3]; // [
																		 // [indexOfParent1, indexOfParent2,indexOfF1_1],
																		 // [indexOfParent,indexOfF1_2],...]
		for (int i = 0; i < outCrossTaxaIndex.size(); i++) {
			for (int j = 1; j < taxaInfo.size(); j++) {
				if (currGenos.taxaName(outCrossTaxaIndex.get(i)).contains(taxaInfo.get(j).get(0))) {
						String p1 = taxaInfo.get(j).get(4);
						String p2 = taxaInfo.get(j).get(6);
						
						p1name = getParentName(p1);
						p2name = getParentName(p2);
						
						if(p1name == null) {
							p1name = "MISSING";
						}if(p2name == null) {
							p2name = "MISSING";
						}						
						
						int p1index = getGenotypeIndex(p1name);
						int p2index = getGenotypeIndex(p2name);
						int childIndex = getGenotypeIndex(currGenos.taxaName(outCrossTaxaIndex.get(i)));
						
						progenyParentMap[i][0] = p1index;
						progenyParentMap[i][1] = p2index;
						progenyParentMap[i][2] = childIndex;
//						System.out.println(childIndex + genos.taxaName(childIndex) + ": " + p1name + " " + p2name);

//						System.out.println(genos.taxaName(p1index) + " " + genos.taxaName(p2index));
				}
			}
		}

		return progenyParentMap;
	}
		
	private static String getParentName(String parentName) {
		for(int i = 0; i < taxaInfo.size(); i++) {
			if(parentName.trim().equals(taxaInfo.get(i).get(1).trim())) {
				return taxaInfo.get(i).get(0);
			}
		}
		
		
		return null;
	}
	
	private static int getGenotypeIndex(String taxaName) {
		int index = 0;
		for(int i = 0; i < genos.numberOfTaxa(); i++) {
			if(genos.taxaName(i).trim().equals(taxaName.trim())) return i;
		}
		
		return index;
	}

	public static GenotypeTable getParentsGenoTable() {
		GenotypeTable pg = getParentsGenotypeTable(genos, taxaInfo);

		return pg;
	}

	private static GenotypeTable getParentsGenotypeTable(GenotypeTable genos, ArrayList<ArrayList<String>> taxaInfo) {
		ArrayList<String> parentNames = new ArrayList<String>();
		ArrayList<String> parentDNANames = new ArrayList<String>();
		TaxaList currTaxa = genos.taxa();
		MutableTaxaList keepParentTaxa = new MutableTaxaList(); // has the current taxa to be fed to
																// GetPollinationType.java (1 sample group per run)

		for (int i = 1; i < taxaInfo.size(); i++) {
			if (!parentNames.contains(taxaInfo.get(i).get(1))) {
				parentNames.add(taxaInfo.get(i).get(1));
				parentDNANames.add(taxaInfo.get(i).get(0));

			}
		}
		parentNames.clear();
		// get parent genotypes from currTaxa

		for (int j = 0; j < parentDNANames.size(); j++) {
			for (int k = 0; k < currTaxa.size(); k++) {
				if (parentDNANames.get(j).matches(currTaxa.get(k).getName())) {
					keepParentTaxa.add(currTaxa.get(k));
				}
			}
		}

		return FilterGenotypeTable.getInstance(genos, keepParentTaxa, false);
	}

	private static GenotypeTable passArrayList(GenotypeTable genos, ArrayList<ArrayList<String>> taxaInfo) {
		GenotypeTable parentGenos = getParentsGenotypeTable(genos, taxaInfo);
		return parentGenos;
	}

	public static ArrayList<GenotypeTable> getFilteredGenotypeTable(String hmpFile, String pedFile) {
		String genotype_file = hmpFile;
		String genotype_file_info = pedFile;

		return readGenotypeTable(genotype_file, genotype_file_info);
	}

	private static ArrayList<GenotypeTable> readGenotypeTable(String genotype_file,
			String genotype_file_info) {
		// READ GENOTYPE FILE INFORMATION // SAMPLES
		taxaInfo = new ArrayList<ArrayList<String>>();
		String line = null;
		int taxaSize = 0; // counter of every ArrayList initialized inside taxa ArrayList
		try {
			BufferedReader buf = new BufferedReader(new FileReader(genotype_file_info));

			while ((line = buf.readLine()) != null) {
				taxaInfo.add(new ArrayList<String>());
				String[] holder = line.split("\t");
				for (int i = 0; i < holder.length; i++) {
					if (holder[i].trim().length() > 0) {
						taxaInfo.get(taxaSize).add(holder[i]);
					} else {
						taxaInfo.get(taxaSize).add(""); // added empty string to maintain nxn matrix
					}
				}
				taxaSize++; // maintain 0-indexing; place counter after adding element inside
							// taxa.get(arraylist)
			}

			buf.close();
		} catch (FileNotFoundException e) {	
			System.out.println("File not found!");
		} catch (IOException e) {
			System.out.println("Error reading file.");
		}

//        for(int i = 0; i < taxaInfo.size(); i++) {
//        		System.out.println(taxaInfo.get(i));
//        }

		// READ GENOTYPE TABLE
		genos = ImportUtils.readFromHapmap(genotype_file);

//		passArrayList(genos, taxaInfo);
		// GET GROUPS
		taxaPerGroup = sampleByGroup(genos, taxaInfo);

		ArrayList<GenotypeTable> listOfGenos = compareGenotypeFiles(genos, taxaInfo, taxaPerGroup);
		return listOfGenos;
	}

	private static HashMap<String, ArrayList<Integer>> sampleByGroup(GenotypeTable genos,
			ArrayList<ArrayList<String>> taxaInfo) {
		// TODO Auto-generated method stub
		TaxaList currTaxa = genos.taxa();
		ArrayList<String> taxaGroups = new ArrayList<String>();
		HashMap<String, ArrayList<Integer>> taxaPerGroup = new HashMap<String, ArrayList<Integer>>();

		for (int i = 0; i < taxaInfo.size(); i++) {
			if (taxaInfo.get(i).size() > 8) {
				if (!taxaGroups.contains(taxaInfo.get(i).get(8))) {
					taxaGroups.add(taxaInfo.get(i).get(8));
				}
			}
		}

		for (int k = 1; k < taxaGroups.size(); k++) {
			String gname = taxaGroups.get(k).toLowerCase();
			ArrayList<Integer> indices = new ArrayList<Integer>();
			for (int l = 1; l < currTaxa.size(); l++) {
				String dname = currTaxa.get(l).getName().toLowerCase();
				if (dname.trim().contains(gname.trim()+"-") && !(taxaPerGroup.containsKey(gname))) {
					indices.add(l);
				}
			}
			taxaPerGroup.put(gname, indices);
		}

//		print
//		for (String name: taxaPerGroup.keySet()){
//           System.out.println(taxaPerGroup.get(name));  
//		} 

		return taxaPerGroup;
	}

	private static ArrayList<GenotypeTable> compareGenotypeFiles(GenotypeTable genos,
			ArrayList<ArrayList<String>> taxaInfo, HashMap<String, ArrayList<Integer>> taxaPerGroup) {
		ArrayList<GenotypeTable> listOfGenos = new ArrayList<GenotypeTable>();
		TaxaList currTaxa = genos.taxa();
		boolean isF1_or_BC1F1 = false;

		for (String groupName : taxaPerGroup.keySet()) {
			ArrayList<Integer> arr = taxaPerGroup.get(groupName);
			MutableTaxaList keepTaxa = new MutableTaxaList(); // has the current taxa to be fed to
																// GetPollinationType.java (1 sample group per run)
			for (int taxonIndex = 0; taxonIndex < arr.size(); taxonIndex++) {
				int t = arr.get(taxonIndex); //////////////////////////////
				isF1_or_BC1F1 = matchGermplasmType(PedigreeFileInfo.genos, taxaInfo, taxonIndex);
				if (isF1_or_BC1F1) {
//					System.out.println("trueee");
					keepTaxa.add(currTaxa.get(t));
					continue;
				}
				keepTaxa.add(currTaxa.get(t));
			}
			GenotypeTable filtered = FilterGenotypeTable.getInstance(genos, keepTaxa, false);
			listOfGenos.add(filtered);
		}
//			

		return listOfGenos;
	}

	private static boolean matchGermplasmType(GenotypeTable genos, ArrayList<ArrayList<String>> taxaInfo,
			int taxonIndex) {
		TaxaList currTaxa = genos.taxa();
		int i = taxonIndex;
		int k = 0;

		while (k < taxaInfo.size()) {
			if (taxaInfo.get(k).get(0).matches(currTaxa.get(taxonIndex).getName())) {
				if (taxaInfo.get(k).get(3).matches("F1") || taxaInfo.get(k).get(3).matches("BC1F1")) {
					return true;
				}
				else return false;
			}
			k++;
		}
		return false;
	}
	
	
	public static GenotypeTable getGenotypeTable() {
		return genos;
	}
	
	public static ArrayList<String> getKeys() {
		ArrayList<String> keys = new ArrayList<String>();
		for(String key: taxaPerGroup.keySet()) {
			keys.add(key);
		}
		return keys;
	}

	public static void main(String[] args) {
//		ArrayList<GenotypeTable> filteredGenos = getFilteredGenotypeTable();
//		getParentsGenoTable();
//		long startTime = System.nanoTime();
//		long endTime = System.nanoTime();
//		System.out.println("Took "+(endTime - startTime) + " ns"); 

//		for(int i = 0; i < filteredGenos.numberOfTaxa(); i++) {
//			System.out.println(filteredGenos.taxaName(i));
//		}
	}

}
