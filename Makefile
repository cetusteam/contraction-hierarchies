include compiler.make

SRC = main.cpp
DEP = $(SRC:%.cpp=%.d) 
main: main.o
	$(CXX) $(LINK) -o $@ $^ $(LIBS)

clean:
	rm main main.o

%.d: %.cpp
	$(SHELL) -ec '$(CXX) -MM $(CXXFLAGS) $< \
                      | sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@; \
                      [ -s $@ ] || rm -f $@'

include $(DEP)
