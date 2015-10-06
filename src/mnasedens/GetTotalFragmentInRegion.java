package mnasedens;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import jbfunctions.BAMInput;
import jbfunctions.TF;

public class GetTotalFragmentInRegion {

	public static void main(String[] args) throws IOException {
		
		// Set the input parameters
		String bam_file_name = args[0];
		String input_feature_file = args[1];
		String output_file_name = args[2];
		int left_win = Integer.parseInt(args[3]);
		int right_win = Integer.parseInt(args[4]);
		int frag_low = Integer.parseInt(args[5]);
		int frag_high = Integer.parseInt(args[6]);
		
		///////////////////////////////////////////////
		// Get the TF HashMap
		ArrayList<TF> tf_list = TF.read_in_tf_list(input_feature_file);
		
		// Open the output buffer
		BufferedWriter output = new BufferedWriter(new FileWriter(output_file_name));
		
		// Write the header
		output.write("name,chr,pos,strand,number_fragments\n");
		
		// Open the bam file
		SAMFileReader bam_file = new SAMFileReader(new File(bam_file_name), new File(bam_file_name + ".bai"));
		
		// Get the total number of reads
		String[] chr_names = new String[16];
		for(int c = 0; c < 16; c++){
			chr_names[c] = Integer.toString(c+1);
		}
		int total_number_of_reads = BAMInput.get_number_aligned_reads(bam_file_name, chr_names);
		
		// Set the DecimalFormat
		DecimalFormat df = new DecimalFormat("#.####");
		
		// Iterate through each chromosome
		for(TF feature : tf_list){
			
			// Get the positions of interest
			int l_win, r_win;
			if(feature.getStrand() == '+'){
				l_win = feature.getPos() - left_win;
				r_win = feature.getPos() + right_win;
			}else{
				l_win = feature.getPos() - right_win;
				r_win = feature.getPos() + left_win;
			}

			// Get the read iterator over this region
			SAMRecordIterator bam_itr = bam_file.queryOverlapping(feature.getChr(), l_win - 50, r_win);
			
			// Set the total number of fragments
			double total_fragments = 0;
			
			// Iterate through each bam_itr
			while(bam_itr.hasNext()){
				
				// Get the read
				SAMRecord read = bam_itr.next();
				
				// Get the read width and midpoint position
				int read_width = read.getInferredInsertSize();
				int midpoint_position = read.getAlignmentStart() + (read_width/2);
				
				// If the read midpoint position is within l_win and r_win and 
				// within the frag_low and frag_high, include in the final count
				if((midpoint_position >= l_win) && (midpoint_position <= r_win) &&
				   (read_width >= frag_low) && (read_width <= frag_high)
				  ){
					total_fragments++;
				}
					
			}	
		

			// Close the bam_itr
			bam_itr.close();
			
			// Normalize to 100E6 reads
			total_fragments *= ((100E6) / (double) total_number_of_reads);
			
			// Write the output
			output.write(feature.toString() + "," + df.format(total_fragments) + "\n");
		
			
		}

		// Close the output buffer
		output.close();
		
		// Close the bam_file buffer
		bam_file.close();
		
	}

}
