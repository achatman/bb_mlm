all: bbts load_data output
	$$(root-config --cxx --cflags --libs) $$VEGAS/lib/libVEGAScommon.so $$VEGAS/lib/libVEGAScoord.so -I$$VEGAS/include/resultsExtractor/include/ -I$$VEGAS/include/common/include/ -I$$VEGAS/include/showerReconstruction2/include -I$$VEGAS/include/coordinates/include/ -I$$VEGAS/include/eventSelection/include/ -D_OAWG -lMinuit -lTreePlayer -Wall bbts.o output.o load_data.o -o bbts

clean:
	rm -rf Bin_* RawData* fitstats* Errors_* summary.csv HIST* Cuts_* Likelihood* Bidirectional* MSWMSL*

clobber: clean
	rm -rf slurm* *.o bbts cerr.out

load_data: 
	$$(root-config --cxx --cflags --libs) -I$$VEGAS/include/resultsExtractor/include/ -I$$VEGAS/include/common/include/ -I$$VEGAS/include/showerReconstruction2/include -I$$VEGAS/include/coordinates/include/ -I$$VEGAS/include/eventSelection/include/ -D_OAWG -lMinuit -lTreePlayer -Wall load_data.cpp -c

output: 
	$$(root-config --cxx --cflags --libs) -I$$VEGAS/include/resultsExtractor/include/ -I$$VEGAS/include/common/include/ -I$$VEGAS/include/showerReconstruction2/include -I$$VEGAS/include/coordinates/include/ -I$$VEGAS/include/eventSelection/include/ -D_OAWG -lMinuit -lTreePlayer -Wall output.cpp -c

bbts: 
	$$(root-config --cxx --cflags --libs) -I$$VEGAS/include/resultsExtractor/include/ -I$$VEGAS/include/common/include/ -I$$VEGAS/include/showerReconstruction2/include -I$$VEGAS/include/coordinates/include/ -I$$VEGAS/include/eventSelection/include/ -D_OAWG -lMinuit -lTreePlayer -Wall bbts.cpp -c
