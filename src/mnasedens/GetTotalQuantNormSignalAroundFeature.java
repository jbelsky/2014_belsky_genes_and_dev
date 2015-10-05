package mnasedens;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import functions.BAMInput;
import functions.TF;


public class GetTotalQuantNormSignalAroundFeature {

	public static double[][] get_chr_density(String[] file_names, String footer2, String ref_bam_file, String chr) throws IOException{
		
		// Get the chr length
		int chr_length = BAMInput.get_chr_length(ref_bam_file, chr);
		
		// Create the storage double
		double[][] mat = new double[file_names.length][chr_length + 1];
		
		// Iterate through each file name
		for(int f = 0; f < file_names.length; f++){
			
			// Read in the relevant file
			BufferedReader input = new BufferedReader(new FileReader(file_names[f] + "_" + chr + "_" + footer2));
			
			// Read in the header
			input.readLine();
			
			// Enter into the storage array
			for(int i = 1; i < mat[0].length; i++){
				String[] line_arr = input.readLine().split(",");
				mat[f][i] = Double.parseDouble(line_arr[2]);
			}
			
			// Close the buffer
			input.close();
			
		}
		
		// Return the storage matrix
		return(mat);
		
		
	}
	
	public static double[] get_total_signal_per_feature(TF t, int l_win, int r_win, double[][] mat){
		
		// Set up the storage output
		double[] total_sig = new double[mat.length];
		
		// Get the position
		int pos = t.getPos();
		
		// Set the start and end positions
		int start = pos - l_win;
		int end = pos + r_win;
		
		// If the strand is negative, flip the start and end
		if(t.getStrand() == '-'){
			start = pos - r_win;
			end = pos + l_win;
		}
		
		// Get the total signal over the range for each dataset
		for(int a = 0; a < total_sig.length; a++){
			for(int p = start; p <= end; p++){
				if(p >= 0 && p < mat[a].length){
					total_sig[a] += mat[a][p];
				}
			}
			
			// Ensure that the total_sig is at least 0
			if(total_sig[a] < 0){
				total_sig[a] = 0;
			}
			
		}
		
		// Return the storage output
		return(total_sig);
		
	}
	
	public static String get_output_string(TF t, double[] total_sig){
		
		// Set the output DecimalFormat
		DecimalFormat df = new DecimalFormat("#.####");
		
		// Set the header string
		String output_string = t.getName() + "," +
							   t.getChr() + "," +
							   t.getPos() + "," +
							   t.getStrand() + ",";
		
		// Iterate through the total_sig
		for(int a = 0; a < total_sig.length; a++){
			String sep = ",";
			if(a == total_sig.length - 1){
				sep = "\n";
			}
			output_string += df.format(total_sig[a]) + sep;
		}
		
		// Return the output string
		return(output_string);
		
	}
	
	public static void main(String[] args) throws IOException {

		// Set the string headers
		String[] dataset_header = {"orc_chip_seq_dm265",
					   "orc_chip_seq_dm287",
					   "input_chip_seq_dm271_orc",
					   "input_chip_seq_dm272_orc"
					  };
				
		// Set the directory
		String dataset_dir = "/data/data2/jab112/2014_mnase_manuscript/datasets/chip_data/orc/quant_normalize/";
		
		// Set the footers
		String footer1 = "density_signal_chr";
		String footer2 = "quant_normalize.csv";
		
		// Set the input feature file name
		String input_feature_file = 
			"/data/data2/jab112/2014_mnase_manuscript/datasets/jab112_yeast_feature_files/replication_origins/" +
			"oridb_acs_feature_file_curated_798_sites_timing_whitehouse_raw_oem.csv";
			
		// Set the window
		int left_win = 400;
		int right_win = 400;
		
		// Set the output file name
		String output_file_name = 
			"/data/data2/jab112/2014_mnase_manuscript/figures/figure_3/figure_3_datasets/" +
			"oridb_acs_feature_file_curated_798_sites_orc_chip_seq_" + 
			"left_win_" + left_win + "bp_right_win_" + right_win + "bp.csv";				  
		
		// Set a reference bam file name
		String bam_file_name = "/data/illumina_pipeline/aligned_experiments/DM242/dm242.bam";
		
		////////////////////////////////////////////////////////////////////////////////////
		// Set the file names
		String[] input_file_name_arr = new String[dataset_header.length];
		for(int h = 0; h < input_file_name_arr.length; h++){
			input_file_name_arr[h] = dataset_dir + dataset_header[h] + "_" + footer1;
		}
		
		// Get the TF HashMap
		HashMap<String, ArrayList<TF>> tf_map = TF.read_in_tf_map(input_feature_file);
		
		// Open the output buffer
		BufferedWriter output = new BufferedWriter(new FileWriter(output_file_name));
		
		// Write the header
		output.write("name,chr,pos,strand,");
		
		for(int h = 0; h < dataset_header.length; h++){
			String sep = ",";
			if(h == dataset_header.length - 1){
				sep = "\n";
			}
			output.write(dataset_header[h] + sep);
		}
		
		// Iterate through each TF
		for(int c = 1; c <= 16; c++){
			
			// Get the double array
			System.out.println("Reading in the data for chr " + c + "...");
			double[][] input_data = get_chr_density(input_file_name_arr, footer2, 
													bam_file_name, Integer.toString(c)
												   );
					
			// Get the tf_list on the chr
			ArrayList<TF> tf_list = tf_map.get(Integer.toString(c));
			
			// Iterate through each feature
			if(tf_list != null){
				
				for(int f = 0; f < tf_list.size(); f++){
				
					// Get the total signal output
					double[] total_signal = get_total_signal_per_feature(tf_list.get(f), left_win, right_win, input_data);
					
					// Write the output
					output.write(get_output_string(tf_list.get(f), total_signal));
				
				}
					
			}
			
		}

		// Close the output buffer
		output.close();
		
	}
	
}
