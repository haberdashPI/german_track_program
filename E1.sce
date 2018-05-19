response_matching = simple_matching;
response_logging = log_active;
active_buttons = 5; #space,y,n,enter,r
default_font_size = 30;
write_codes = true; # send codes to biosemi box (via the parallel port)
stimulus_properties = block_type,string,block_idx,number,trial_order,number,trial_file,number;
event_code_delimiter = ";";

begin;

wavefile { filename = ""; preload = false; } sound_file;
text { 
	caption = "temp"; 
	preload = false;
	font_size = 30; 
	max_text_width = 1350;
} message;

trial {
   nothing { default_port_code = 64; };
	duration = 100;
} start_recording;

trial {
	nothing { default_port_code = 128; };
	duration = 100;
} stop_recording;

trial {    
	stimulus_event {
		sound { 
			wavefile sound_file;
		};   
		port_code = 16;
		port = 1;
		code = "guided_trial";
	};	
	LOOP $i 1;
	stimulus_event {
		picture { 
			text message; 
			x = 0; y = 0;
		} guide_msg;
	};	
	ENDLOOP;
} guide_trial;

trial {    
	stimulus_event {
		sound {
			wavefile sound_file;
		};   
		time = 0;
		code = "train;0;0;0";
		port_code = 16;
		port = 1;
	} se;	
} sound_trial;

trial {
	trial_type = specific_response;
	terminator_button = 1;
	trial_duration = forever;
	picture { 
		text message; 
		x = 0; y = 0;
	};
} message_trial;

trial {
	trial_type = specific_response;
	terminator_button = 2,3;
	trial_duration = forever;
	picture { 
		text message; 
		x = 0; y = 0;
	};
} question_trial;

trial {
	trial_type = specific_response;
	terminator_button = 1,5;
	trial_duration = forever;
	picture { 
		text message; 
		x = 0; y = 0;
	};
} answer_trial;

trial {
   all_responses = false;
	trial_type = specific_response;
	terminator_button = 1;
	trial_duration = forever;

   picture { text { caption = "Rest for a minute..."; }; x = 0; y = 0; };
	time = 0;
   duration = 60000;
	stimulus_time_in = 1000;
	
   picture { text message; x = 0; y = 0; };
	time = 60100;
   target_button = 1;
	stimulus_time_in = 0;
   stimulus_time_out = never;
} rest_trial;

trial {
	picture {} default;
} blank_trial;

begin_pcl;

#--------------------------------------------------------------------------
# EXPERIMENT CONFIGURATION
#--------------------------------------------------------------------------
string read_path = "stim_wav\\";
string train_path = "stim_training\\";
string record_path = "stim_records\\";
string save_path = "sbj_logs\\";
string term_text = "";
string term_text_this;

# experiment params
#int group = random(0,1); #0 does direction first, 1 does speaker first
int group = 1;
bool do_training = true; # TEMPORARY, normally, true;
bool do_training_loud = true; # TEMPORARY, normally, true;
bool do_feat_orientation = true;
bool do_obj_orientation = true;
string current_path;
string current_path_txt;
array<string> messages[0];
array<string> train_messages[0];
array<double> timestamps[0]; 
array<double> is_control[0]; 
array<double> speaker[0]; 
array<double> direction[0]; 
#-----------------------------------mL fR mR fL fR mL
#array<int> train_ordering_loud[6] = {2,10,17,22,34,1};
array<int> train_ordering_loud[4] = {2,10,17,22};
#-------------------------------mLfR C  mR C  C  fL mL fR C
array<int> train_ordering[10] = {4,12,13,19,28,37,24,32,36,40};
array<int> orient1[2] = {23,31};
array<int> orient2[2] = {20,33};
array<int> ordering[0];
int num_trial_per_block = 50; # TEMPORARY normally 50
int until_trial_idx = 50; # TEMPORARY normally 50
int break_trial_idx = 17; # TEMPORARY normally 17
bool skip_training = false;
int num_sal_trials;
int trial_answer;
bool is_correct_answer = false;

# performance params
stimulus_data last_stim;
response_data curr_resp;
int num_resp;
double num_correct;
double num_false_correct;
double num_hits;
double num_hits_block;

int msg_train_start_1 = 1;
int msg_train_start_2 = 2;
int msg_train_start_3 = 3;
int msg_train_start_4 = 14;
int msg_train_middle_1 = 15;

int msg_test_start_1 = 4;

int msg_train_right_1 = 5;
int msg_train_right_2 = 6;
int msg_test_right_1 = 7;
int msg_test_right_2 = 8;

int msg_train_male_1 = 9;
int msg_train_male_2 = 10;
int msg_test_male_1 = 11;
int msg_test_male_2 = 12;
int msg_done = 13;

include "E1_subs.pcl";
read_messages("messages.txt");
read_train_messages("train_messages.txt");

#--------------------------------------------------------------------------
# RUN EXPERIMENT
#--------------------------------------------------------------------------

print_termfile_line("Things are working, I am starting the experiment.");
print_termfile_line("Writing to: " + logfile_directory);
print_termfile_line("Subject id: " + logfile.subject());
print_termfile_line("Subject group: " + string(group));
print_termfile_line(date_time());
message_showing(messages[msg_train_start_1]);
message_showing(messages[msg_train_start_2]);
message_showing(messages[msg_train_start_3]);
message_showing(messages[msg_train_start_4]);
print_termfile_line("I am starting training.");
	
string block_type = "train";
int until_block_idx = 3; 
int block_count = 0; 
int file_idx = -1;
loop int block_idx = 1 until block_idx > until_block_idx
begin
	
	term_text_this = "";
	print_termfile_line(date_time());
	print_termfile_line("Doing block: " + string(block_idx) + ". Block count: " + string(block_count));
	num_correct = 0;
	num_false_correct = 0;
	num_hits_block = 0;
	num_sal_trials = 0;
	
	# get new random order for trials - this block only
	ordering.resize(0); 
	int rep_num = int(ceil(num_trial_per_block/10));
	loop int j = 1 until j > int(ceil(until_trial_idx/10))
	begin
		loop int i = j until i > num_trial_per_block
		begin
			ordering.add(i);
			i = i + rep_num;
		end;
		j = j + 1;
	end;
	ordering.shuffle();
	
	# present block trials --------------------------------------------------
	loop int trial_idx = 1 until trial_idx > until_trial_idx
	begin
		# trigger eeg recording
		start_recording.present();
			
		if do_training_loud then
			file_idx = train_ordering_loud[trial_idx];
			current_path = train_path + "trial_" + string(file_idx);
			read_timestamps("stim_training/target_info_all.txt");
		elseif do_training then
			file_idx = train_ordering[trial_idx];
			current_path = train_path + "trial_" + string(file_idx);
			read_timestamps("stim_training/target_info_all.txt");
		elseif (block_idx==3 && do_obj_orientation) || (block_idx==2 && do_feat_orientation) then
			file_idx = trial_idx;
			current_path = train_path + "trial_" + string(file_idx);
		else
			file_idx = ordering[trial_idx];
			current_path = read_path + "trial_" + string(file_idx);
			if block_idx==1 then
				read_timestamps("target_info_all.txt");
			elseif block_idx==3 then
				read_timestamps("target_info_obj.txt");
			else
				read_timestamps("target_info_dir.txt");
			end;
		end;
		print_termfile("Trial: " + string(trial_idx) + " stim " + string(file_idx) + ".\t");
		
		# trial presentation -------------------------------------------------
		soundtrial_button_activation();
		sound_file.set_filename(current_path + ".wav");
		sound_file.load();
		se.set_event_code(block_type+";"+string(block_idx)+";"+string(trial_idx)+";"+string(file_idx));
		if block_type == "train" then
			se.set_port_code(16);
		else
			se.set_port_code(32);
		end;
		sound_trial.present(); # PRESENT IS HERE
		
		# for orientation in cue block
		if block_idx==3 && do_obj_orientation then
			if trial_idx==1 then
				message_showing("Press '4' to hear another example trial.");
				trial_idx = trial_idx + 1;
				continue;
			else
				block_type = "object";
				do_obj_orientation = false;
				message_showing(messages[msg_test_male_1]);
				message_showing(messages[msg_test_male_2]);
				term.print_line("");
				trial_idx = 1;
				continue;
			end;
		elseif block_idx==2 && do_feat_orientation then
			if trial_idx==1 then
				message_showing("Press '4' to hear another example trial.");
				trial_idx = trial_idx + 1;
				continue;
			else
				block_type = "feature";
				do_feat_orientation = false;
				message_showing(messages[msg_test_right_1]);
				message_showing(messages[msg_test_right_2]);				
				term.print_line("");
				trial_idx = 1;
				continue;
			end;
		end;
		
		# ask the saliency q and check response ------------------------------
		questiontrial_button_activation();
		message.set_caption(msg_q + "\n\n");
		message.redraw();
		question_trial.present();
		# check that they give a valid answer before moving on 
		loop int last_bttn = response_manager.last_response() until last_bttn==4
		begin
			trial_answer = response_manager.last_response();
			if (trial_answer == 2) then
				message.set_caption(msg_q + "yes\n\n" + msg_q_moveon);
			else
				message.set_caption(msg_q + "no\n\n" + msg_q_moveon);
			end;
			message.redraw();
			question_trial.set_terminator_buttons({2,3,4});
			question_trial.present();
			last_bttn = response_manager.last_response();
		end;
		question_trial.set_terminator_buttons({2,3});
		blank_trial.present();
		# then check if that answer is correct
		if trial_answer==2 && is_control[file_idx]==0 then
			num_correct = num_correct + 1;
		elseif trial_answer==2 && is_control[file_idx]==1 then
			num_false_correct = num_false_correct + 1;
		end;
		if is_control[file_idx]==0 then
			num_sal_trials = num_sal_trials + 1;
		end;
		
		# wrapping up the trial ----------------------------------------------
		answertrial_button_activation();
		if do_training || do_training_loud then
			print_termfile_line("Trial "+string(1-is_control[file_idx])+" answer "+string(3-trial_answer)
				+" speaker "+string(speaker[file_idx])+" direction "+string(direction[file_idx]));
			string this_text;
			is_correct_answer = false;
			if (trial_answer==2 && is_control[file_idx]==0) then
				this_text = "Correct!\nThere WAS a different pitch segment in this clip.";
				is_correct_answer = true;
			elseif (trial_answer==2 && is_control[file_idx]==1) then
				this_text = "Wrong!\nThere WAS NO different pitch segment in this clip.";
			elseif (trial_answer==3 && is_control[file_idx]==1) then
				this_text = "Correct!\nThere WAS NO different pitch segment in this clip.";
				is_correct_answer = true;
			elseif (trial_answer==3 && is_control[file_idx]==0) then
				this_text = "Wrong!\nThere WAS a different pitch segment in this clip.";
			end;
			string this_text1;
			string this_text2;
			if is_control[file_idx]==0 then
				this_text = this_text + "\n\n" + train_messages[file_idx];
				this_text1 = this_text + "\n\nPress '4' to hear the trial again, with the pitch segment location shown to you.";
				this_text2 = this_text + "\n\nPress '3' to listen to the clip again with the pitch segment location shown to you.";
			else
				this_text = this_text + "\n\n";
				this_text1 = this_text + "\n\nPress '4' to hear the trial again.";
				this_text2 = this_text + "\n\nPress '3' to listen to the clip again.";
			end;
			this_text2 = this_text2 + "\n\nPress '4' to start the next trial.";
			if do_training_loud then
				message.set_caption(this_text1);
			else
				message.set_caption(this_text2);
			end;
			message.redraw();
			answer_trial.present();
			# check whether they want to replay the trial 
			bool did_once = false;
			loop int last_bttn = response_manager.last_response() 
			until last_bttn==1 && (!do_training_loud || (do_training_loud&&did_once))
			begin
				print_termfile_line("repeating trial");
				disable_buttons();
				blank_trial.present();
				message.set_caption("Here");
				message.redraw();
				if direction[file_idx]==3 then
					guide_msg.set_part_x(1,500);
				else
					guide_msg.set_part_x(1,-500);
				end;
				if is_control[train_ordering[trial_idx]]==1 then
					sound_trial.present();
				else
					guide_trial.get_stimulus_event(2).set_time(int(timestamps[file_idx]*1000)); #in msec
					guide_trial.get_stimulus_event(2).set_duration(1000);
					guide_trial.present();
				end;
				answertrial_button_activation();
				message.set_caption("Your last answer was: " + this_text2);
				message.redraw();
				answer_trial.present();
				last_bttn = response_manager.last_response();
				did_once = true;
			end;
			blank_trial.present();
			if !is_correct_answer then
				term.print_line("Adding this trial to the end of the queue.");
				if do_training_loud then
					train_ordering_loud.add(train_ordering_loud[1]);
				else
					train_ordering.add(train_ordering[1]);
				end;				
			end;
			if do_training_loud then
				remove(train_ordering_loud,1);
			else
				remove(train_ordering,1);
			end;				
			if do_training_loud && train_ordering_loud.count()==0 then
				message_showing(messages[msg_train_middle_1]);
				do_training_loud = false;
			elseif do_training && train_ordering.count()==0 then
				message_showing(messages[msg_test_start_1]);
				do_training = false;
				block_type = "test";
				num_correct = 0;
				num_false_correct = 0;
				num_hits_block = 0;
				num_sal_trials = 0;
				print_termfile_line("Starting experiment section.");
				print_termfile_line("Doing block: " + string(block_idx));
			end;
		else
			print_termfile_line("Trial "+string(1-is_control[file_idx])+" answer "+string(3-trial_answer)
				+" speaker "+string(speaker[file_idx])+" direction "+string(direction[file_idx]));
				
			if mod(trial_idx,break_trial_idx)==0 && trial_idx!=until_trial_idx then
				print_termfile("Giving a break with the rest trial.\n");
				message.set_caption("Press '4' to start the next trial." + rest_msgs[block_idx]);
				message.redraw();
				soundtrial_button_activation();
				stop_recording.present();
				rest_trial.present();
				blank_trial.present();
				start_recording.present();
			end;
			
			trial_idx = trial_idx + 1;
		end;
		
	end;
	
	# calculate block stats for sanity --------------------------------------
	double hr; 
	double fr;
	if num_sal_trials==0 then
		hr = 0;
	else
		hr = round(num_correct/num_sal_trials,2);
	end;
	if until_trial_idx==num_sal_trials then
		fr = 0;
	else
		fr = round(num_false_correct/(until_trial_idx-num_sal_trials),2);
	end;
	print_termfile_line("Block analytics: -------------------------------");
	print_termfile_line("Hit rate: "+string(hr));
	print_termfile_line("False rate: "+string(fr));
	print_termfile_line("------------------------------------------------");
  
	# wrapping up the block -------------------------------------------------
	block_count = block_count+1;
	output_file out = new output_file;
	out.open(save_path + "terminal_" + logfile.subject() + ".txt"); 
	out.print(term_text);  
	out.close();
	if block_idx==1 then
		if group==0 then
			block_idx=2;
		else
			block_idx=3;
		end;
	else
		if group==0 then
			block_idx=3;
		else
			block_idx=2;
		end;
	end;
	string this_msg;
	if block_idx==3 then
		block_type = "object";
		this_msg = messages[msg_train_male_1];
	elseif block_idx==2 then 
		block_type = "feature";
		this_msg = messages[msg_train_right_1];
	end;
	if block_count==1 then
		message_showing(msg_block2+this_msg);
	elseif block_count==2 then
		message_showing(msg_block3+this_msg);
	else
		break;
	end;
	
	if block_idx==3 then
		message_showing(messages[msg_train_male_2]);
	elseif block_idx==2 then
		message_showing(messages[msg_train_right_2]);
	end;
end;

stop_recording.present();

print_termfile_line("\nAll done!");
print_termfile_line(date_time());
message_showing(messages[msg_done]);
output_file out = new output_file;
out.open(save_path + "terminal_" + logfile.subject() + ".txt"); 
out.print(term_text);  
out.close();









