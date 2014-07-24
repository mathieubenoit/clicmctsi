CC=g++
ROOTCFLAGS=`root-config --cflags`
CFLAGS=-c -Wall $(ROOTCFLAGS)
LDFLAGS=`root-config --evelibs`
SOURCES= AllPixDigitAnimation.cpp Charge.cpp ElectricField.cpp Electron.cpp Hole.cpp Interaction.cpp MCTSi.cpp PixelGeometry.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=MCTSi

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm ./*.o
	rm ./MCTSi