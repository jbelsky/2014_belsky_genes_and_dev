package mnasedens;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import jbfunctions.BAMInput;
import jbfunctions.TF;


public class GetQuantNormSignalAroundFeature {

	public static double[] get_chr_density(String file_name, String ref_bam_file, String chr) throws IOException{
		
		// Get the chr length
		int chr_length = BAMInput.get_chr_length(ref_bam_file, chr);
		
		// Create the storage double
		double[] mat = new double[chr_length + 1];
		
		// Read in the relevant file
		BufferedReader input = new BufferedReader(new FileReader(file_name));
			
		// Read in the header
		input.readLine();
			
		// Enter into the storage array
		for(int i = 1; i < mat.length; i++){
			String[] line_arr = input.readLine().split(",");
			mat[i] = Double.parseDouble(line_arr[2]);
		}
			
		// Close the buffer
		input.close();
		
		// Return the storage matrix
		return(mat);		
		
	}
	
	public static double[] get_feature_output(TF t, double[] mat, int win){
		
		// Set up the storage output
		double[] total_sig = new double[mat.length];
		
		// Set the signal idx
		int sig_idx = 0;
		
		// Get the position
		int pos = t.getPos();
		
		// Set the start and end positions
		int start = pos - win;
		int end = pos + win;
		
		// Iterate through the positions
		for(int p = start; p <= end; p++){
			
			// If p goes outside the boundary of the chromosome ends, enter a -1
			if(p < 1 || p >= total_sig.length){
				total_sig[sig_idx] = -1;
			}else{
				total_sig[sig_idx] = mat[p];
			}
			
			// Increment the sig_idx
			sig_idx++;
			
		}

		
		// Return the storage output
		return(total_sig);
		
	}
		
	public static void main(String[] args) throws IOException {

		// Set the string header
		String dataset_header = args[0];
		String dataset_footer = args[1];		
		String input_feature_file = args[2];
		String output_file_name = args[3];
		String bam_file_name = args[4];
		int win = Integer.parseInt(args[5]);
				
		////////////////////////////////////////////////////////////////////////////////////
		
		// Get the TF HashMap
		HashMap<String, ArrayList<TF>> tf_map = TF.read_in_tf_map(input_feature_file);
		
		// Open the output buffer
		BufferedWriter output = new BufferedWriter(new FileWriter(output_file_name));
		
		// Write the header
		output.write("name,chr,pos,strand,");
		String sep = ",";
		for(int i = -win; i <= win; win++){
			if(i == win){
				sep = "\n";
			}
			output.write(i + sep);
		}
		
		// Iterate through each chromosome
		for(int c = 1; c <= 16; c++){
			
			// Set the filename
			String file_name = dataset_header + c + dataset_footer;
			
			// Get the double array
			System.out.println("Reading in the data for chr " + c + "...");
			double[] input_data = get_chr_density(file_name, bam_file_name, Integer.toString(c));
					
			// Get the tf_list on the chr
			ArrayList<TF> tf_list = tf_map.get(Integer.toString(c));
			
			// Iterate through each feature
			if(tf_list != null){
				
				for(int f = 0; f < tf_list.size(); f++){
				
					// Get the TF Feature
					TF feature = tf_list.get(f);
					
					// Get the total signal output
					double[] feature_signal = get_feature_output(feature, input_data, win);
					
					// Write the output
					feature.write_reads_output(output, feature_signal, 1);
				
				}
					
			}
			
		}

		// Close the output buffer
		output.close();
		
	}
	
}
