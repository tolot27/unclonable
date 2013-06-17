import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;


public class lesenSam {

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
	
	
	
	
public static HashMap einlesenSam(String file_name, String contig_id){
		
		File file = new File(file_name);
		HashMap<String,String> map = new HashMap<String, String>();
			
		try {
			BufferedReader buf = new BufferedReader(new FileReader(file));

			String zeile;
			while ((zeile = buf.readLine()) != null) {
				
				
				
				if (!zeile.startsWith("@")) {
				
					
				String[] zeileArray = zeile.split("\t");
				if(zeileArray[2].equals(contig_id)){
				if(Integer.parseInt(zeileArray[3])!=0){
				int stop = Integer.parseInt(zeileArray[3])+zeileArray[9].length()-2;
				
				if(Integer.parseInt(zeileArray[8])>0){
				map.put(zeileArray[0], "F$"+(Integer.parseInt(zeileArray[3])-1)+"$"+ stop);	
					
				}
				else{
					map.put(zeileArray[0], "R$"+(Integer.parseInt(zeileArray[3])-1)+"$"+ stop);	
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
		return map;
		}





}
