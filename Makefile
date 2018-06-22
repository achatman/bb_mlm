all: bbts load_data output
	$$(root-config --cxx --cflags --libs) $$VEGAS/common/lib/libSP24sharedLite.so $$VEGAS/resultsExtractor/lib/libStage6shared.so $$VEGAS/coordinates/lib/vaCoordShared.so $$VEGAS/showerReconstruction2/lib/libStage4.so -I$$VEGAS/resultsExtractor/include -I$$VEGAS/common/include -I$$VEGAS/showerReconstruction2/include -I$$VEGAS/coordinates/include -I$$VEGAS/eventSelection/include -D_OAWG -lMinuit -lTreePlayer load_data.o bbts.o output.o -o bbts

clean:
	rm -rf Bin_* RawData* fitstats* Errors_* summary.csv HIST* Cuts_* Likelihood* Bidirectional* MSWMSL*

clobber: clean
	rm -rf slurm* *.o bbts cerr.out

load_data:
	$$(root-config --cxx --cflags --libs) -I$$VEGAS/resultsExtractor/include -I$$VEGAS/common/include -I$$VEGAS/showerReconstruction2/include -I$$VEGAS/coordinates/include -I$$VEGAS/eventSelection/include -D_OAWG -lMinuit -lTreePlayer -Wall -Werror load_data.cpp -c

output:
	$$(root-config --cxx --cflags --libs) -I$$VEGAS/resultsExtractor/include -I$$VEGAS/common/include -I$$VEGAS/showerReconstruction2/include -I$$VEGAS/coordinates/include -I$$VEGAS/eventSelection/include -D_OAWG -lMinuit -lTreePlayer -Wall -Werror output.cpp -c

bbts:
	$$(root-config --cxx --cflags --libs) -I$$VEGAS/resultsExtractor/include -I$$VEGAS/common/include -I$$VEGAS/showerReconstruction2/include -I$$VEGAS/coordinates/include -I$$VEGAS/eventSelection/include -D_OAWG -lMinuit -lTreePlayer -Wall -Werror bbts.cpp -c
