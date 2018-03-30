all: minuit

clean:
	rm -rf bbts BIN* raw_data* fitstats* error* summary.csv HIST* cuts* BB_Likelihood* Std_Likelihood*

clobber: clean
	rm -rf slurm*

minuit:
	$$(root-config --cxx --cflags --libs) -I$$VEGAS/resultsExtractor/include -I$$VEGAS/common/include -I$$VEGAS/showerReconstruction2/include -I$$VEGAS/coordinates/include -I$$VEGAS/eventSelection/include $$VEGAS/common/lib/libSP24sharedLite.so $$VEGAS/resultsExtractor/lib/libStage6shared.so $$VEGAS/coordinates/lib/vaCoordShared.so $$VEGAS/showerReconstruction2/lib/libStage4.so -D_OAWG -lMinuit -lTreePlayer -Wall -Werror bbts.cpp -c
	$$(root-config --cxx --cflags --libs) -I$$VEGAS/resultsExtractor/include -I$$VEGAS/common/include -I$$VEGAS/showerReconstruction2/include -I$$VEGAS/coordinates/include -I$$VEGAS/eventSelection/include $$VEGAS/common/lib/libSP24sharedLite.so $$VEGAS/resultsExtractor/lib/libStage6shared.so $$VEGAS/coordinates/lib/vaCoordShared.so $$VEGAS/showerReconstruction2/lib/libStage4.so -D_OAWG -lMinuit -lTreePlayer -Wall -Werror load_data.cpp -c
	g++ bbts.o load_data.o -o bbts

