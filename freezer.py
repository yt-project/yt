from bbfreeze import Freezer
f = Freezer("yt-0.3", includes=("pytz","numpy.oldnumeric.fft","pytz.timezone"))
f.addScript("scripts/reason")
f.addScript("scripts/yt")
f()    # starts the freezing process

