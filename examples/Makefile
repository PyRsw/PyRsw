RM = rm


# Usage: make videos RUN=<run_name>
# 		<run_name> gives the directory i.e. Outputs/<run_name>
# Result: Creates .mp4 videos from the frames.
videos:
ifndef RUN
	$(error RUN is not set 2)
endif	
	$(foreach var, \
		$(subst _00000.png, , \
			$(subst Outputs/$(RUN)/Frames/, , \
				$(wildcard Outputs/$(RUN)/Frames/*_00000.png))), \
		ffmpeg -i Outputs/$(RUN)/Frames/$(var)_\%05d.png -q 1 -y Outputs/$(RUN)/$(var)_animation.mp4;)

# Usage: make clean RUN=<run_name>
# 		<run_name> gives the directory i.e. Outputs/<run_name>
# Result: Deletes files associated with run. 
clean:
ifndef RUN
	$(error RUN is not set 2)
endif	
	$(RM) -rf Outputs/$(RUN)

# Usage:  make clean_all
# Result: Deletes ALL outputs.
clean_all:
	$(RM) -rf Outputs/

help:
	@ echo "The following make targets are available:"
	@ echo "    videos RUN=<run_name>"
	@ echo "        Produces mp4 files of the outputs"
	@ echo "        <run_name> specifies the run"
	@ echo "    clean RUN=<run_name>"
	@ echo "        Deletes the files associated with the run."
	@ echo "        <run_name> specifies the run"
	@ echo "    clean_all "
	@ echo "        Deletes ALL output files."
