## Some questions & suggestions

Hey Kieran, while the 7k shots are running though your script, here are some thoughts I have after parsing your logbook and code:

### Documentation

Your lab book is good, tracking the parameters that you're tuning. 
I recommend also including parameters that you aren't (intentionally) changing. Things like Waveplate angles (at least at the start of a series of experiemnts), if they're set once. Things that one would prefer to set and forget - like laser power - tend not to behave themselves, so logging these carefully is recommended also. The lab book log should contain sufficient information for someone (you!) to reproduce that experiment even after a case of retrograde amnesia. Bonus points for including reasoning behind parameter choices (calculated/measured optima, though process behind method, etc). 

### Method

It's no small task getting an experimental sequence built up - kudos, there's a lot of detail here that you've managed to work through.

Indeed, holding the beam on for an entire sequence gives maximum interrogation time in principle. Some thoughts:
	* How do you know the laser is aligned at all stages of the sequence? Likely isn't, as final trap is shunted somewhat.
	* What is the Zeeman shift at each of these stages? Is the scan wide enough to see them? How does the expected linewidth compare to the splitting between the resonance at different stages?
		* Have you measured the fields at these stages? We know the simulation is only good within a factor of order 1.
	* In your log book you suggest that holding longer is possible; have you tried? What timescale do you expect to reach?
		Why not just hold longer here?
How are you using - or intend to use - the analog import log? Have you ensured the recording is functioning as you want? I can't read any useful information from the diagnostics.



### Analysis

In your analysis, I recommend creating and analysing multiple signals in the same vein as you are currently doing with atom number. If you create a git repository, I'll be happy to follow the development and give you unsolicited advice ;). If you don't want to make it public, that's fine, I think Github has free private repos for students. Otherwise, GitLab definitely has free private repositories. 

Unsurprisingly, the analog import subroutine takes aaaages. Even after the importing, looks like processing is very costly on both RAM and CPU.  UPDATE Almost an hour for the 7000 shot run! And ran out of memory! This is a problem! 
It might be worth writing a little program that post-processes analog imports on the fly during data collection, which will save time in analysis. Alternatively, rework the analog import/preprocess to reduce the RAM load, possibly at the cost of disk access. Consider batch-processing analog import files.  
