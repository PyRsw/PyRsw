videos:
	$(foreach var, \
		$(subst _00000.png, , \
			$(subst Outputs/$(RUN)/Frames/, , \
				$(wildcard Outputs/$(RUN)/Frames/*_00000.png))), \
		ffmpeg -i Outputs/$(RUN)/Frames/$(var)_\%05d.png -q 1 -y Outputs/$(RUN)/$(var)_animation.mp4;)
