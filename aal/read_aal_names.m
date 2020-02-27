function [region_names] = read_aal_names(names_file)
%% Read in the Harvard-Oxford region names

if nargin < 1
    names_file = '/raizadaUsers/mwu34/Documents/MATLAB/spm12/toolbox/aal/atlas/AAL2.xml';
end

all_xml_lines = textread(names_file,'%s','delimiter','\n');

num_lines = length(all_xml_lines);

region_names = struct();
roi_ctr = 0;
for line_num = 1:num_lines,
   this_line = all_xml_lines{line_num};
   
   %%% The lines with regions have the word "index" in them, e.g.
   %%% <label index="17" x="28" y="39" z="63">Superior Parietal Lobule</label>
   starting_char_of_index_string = strfind(this_line,'index');

   if ~isempty(starting_char_of_index_string),  % If we did actually find "index"
      %%% The number after index starts 7 chars later
      number_string = str2num(this_line(starting_char_of_index_string(1) + [6:9]));
      if number_string > 8999
          break;
      end
      roi_ctr = roi_ctr + 1;
      %%% The atlas-image intensity value corresponding to this index is
      %%% one greater than the index number
      
      %%% The name of the region is between triangular parentheses
      opening_triangle_parenth_positions = strfind(this_line,'<name');         
      closing_triangle_parenth_positions = strfind(this_line,'/name>'); 
      
      this_region_name = this_line((opening_triangle_parenth_positions+6): ...
                                   (closing_triangle_parenth_positions-4));
      left_right = this_line(closing_triangle_parenth_positions-2);
      region_names(roi_ctr).name = this_region_name;
      region_names(roi_ctr).index = number_string;
      region_names(roi_ctr).left_right = left_right;
      
   end;
end;