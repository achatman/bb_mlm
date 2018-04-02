all: bbts load_data
	$$(root-config --cxx --cflags --libs) $$VEGAS/common/lib/libSP24sharedLite.so $$VEGAS/resultsExtractor/lib/libStage6shared.so $$VEGAS/coordinates/lib/vaCoordShared.so $$VEGAS/showerReconstruction2/lib/libStage4.so -I$$VEGAS/resultsExtractor/include -I$$VEGAS/common/include -I$$VEGAS/showerReconstruction2/include -I$$VEGAS/coordinates/include -I$$VEGAS/eventSelection/include -D_OAWG -lMinuit -lTreePlayer load_data.o bbts.o -o bbts

clean:
	rm -rf bbts BIN* raw_data* fitstats* error* summary.csv HIST* cuts* BB_Likelihood* Std_Likelihood*

clobber: clean
	rm -rf slurm*

load_data:
	$$(root-config --cxx --cflags --libs) -I$$VEGAS/resultsExtractor/include -I$$VEGAS/common/include -I$$VEGAS/showerReconstruction2/include -I$$VEGAS/coordinates/include -I$$VEGAS/eventSelection/include -D_OAWG -lMinuit -lTreePlayer -Wall -Werror load_data.cpp -c

bbts:
	$$(root-config --cxx --cflags --libs) -I$$VEGAS/resultsExtractor/include -I$$VEGAS/common/include -I$$VEGAS/showerReconstruction2/include -I$$VEGAS/coordinates/include -I$$VEGAS/eventSelection/include -D_OAWG -lMinuit -lTreePlayer -Wall -Werror bbts.cpp -c
