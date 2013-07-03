import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;


public class lesenSam {
	static int zaehler;
	
	public static int getZaehler() {
		return zaehler;
	}
	
	public lesenSam(){
		zaehler = 0;
	}

	public static List<String[]> info_contig(String file_name){
		
		List<String[]> list = new LinkedList<String[]>();
		
		File file = new File(file_name);
		
		try {
			BufferedReader buf = new BufferedReader(new FileReader(file));

			String zeile;
			while ((zeile = buf.readLine()) != null) {
				
				if (zeile.startsWith("@")) {
				 String sn = zeile.split("\t")[1];
					String contig_id = sn.split(":")[1];
					String contig_laenge = zeile.split("\t")[2].split(":")[1];
				
					String[] array = new String[2];
					array[0] = contig_id;
					array[1]= contig_laenge;
					list.add(array);
				}
			}
			buf.close();
			
		}
		catch (Exception e) {
			e.printStackTrace();
			
		}
		return list;
		
		
		
		
		
	}
	
	public static int[] nichtGemappt(String file_name, String contig_id, int contig) throws IOException{
		
		File file = new File(file_name);
		int unmapped =0;
		int mapped = 0;
		
		
		 BufferedWriter Ausgabe_Ungemappt = new BufferedWriter(new FileWriter(new File(
					"nichtGemappt.txt")));
		
						
		
		
			
		try {
			BufferedReader buf = new BufferedReader(new FileReader(file));

			String zeile;
			while ((zeile = buf.readLine()) != null) {
				
				
				
				if (!zeile.startsWith("@")) {
				
					
				String[] zeileArray = zeile.split("\t");
				if(Integer.parseInt(zeileArray[3])==0){
					unmapped++;
					
					
					Ausgabe_Ungemappt.write(zeileArray[0]+"\n");
				
					
				
				}
				else{
					if(zeileArray[2].equals(contig_id)){
					mapped++;
					}
				}
			
				}
			}
			buf.close();
			
		}
		catch (Exception e) {
			e.printStackTrace();
			
		}
		
		Ausgabe_Ungemappt.close();
		
		int[] rueck = new int[2];
		
		rueck[0] = mapped;
		rueck[1] = unmapped;
		
		return rueck;
	}
	
	
public static HashMap einlesenSam(String file_name, String contig_id, int contig, String Organismus_Name) throws IOException{
		
		File file = new File(file_name);
		HashMap<String,String> map = new HashMap<String, String>();
		
		String id = contig_id;
			
		BufferedWriter Ausgabe_SplitRead = new BufferedWriter(new FileWriter(new File(
				"SplitRead_"+Organismus_Name+":"+id+".txt")));
		
		try {
			BufferedReader buf = new BufferedReader(new FileReader(file));

			String zeile;
			while ((zeile = buf.readLine()) != null) {
				
				
				
				if (!zeile.startsWith("@")) {
				
					
				String[] zeileArray = zeile.split("\t");
			
				if(!zeileArray[2].equals("*")){
				if(zeileArray[2].equals(contig_id)){
				if(Integer.parseInt(zeileArray[3])!=0){
				int stop = Integer.parseInt(zeileArray[3])+zeileArray[9].length()-2;
				
				if(!map.containsKey(zeileArray[0])){
				if(Integer.parseInt(zeileArray[8])>0){
				map.put(zeileArray[0], "F$"+(Integer.parseInt(zeileArray[3])-1)+"$"+ stop);	
					
				}
				else{
					map.put(zeileArray[0], "R$"+(Integer.parseInt(zeileArray[3])-1)+"$"+ stop);	
				}
				}else{
					Ausgabe_SplitRead.write(zeileArray[0]+"\n");
					zaehler++;
				}
				}
				}
				}
				}
			}
			buf.close();
			
		}
		catch (Exception e) {
			e.printStackTrace();
			
		}
		
		Ausgabe_SplitRead.close();
		return map;
		}





}
