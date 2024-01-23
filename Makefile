#
# WARNING : g++ and root-system are requiered to compile this example

#
# If pkg-config isn't installed on your system, comment the following lines and
# set the fasterac flags as indicated with your own paths:
#
# FASTERAC_CFLAGS = -I/usr/include
# FASTERAC_LIBS   = -L/usr/lib -lfasterac
#
FASTERAC_CFLAGS = $(shell pkg-config  --cflags libfasterac)
FASTERAC_LIBS   = $(shell pkg-config  --libs   libfasterac)
#

ROOT_CFLAGS     = $(shell root-config --cflags)
ROOT_LIBS       = $(shell root-config --libs)

# Options de liage
LDFLAGS = $(shell root-config --ldflags)

CXX        = g++
CXXFLAGS    = ${FASTERAC_CFLAGS} ${ROOT_CFLAGS} -pthread -Iinclude
LIBS      = ${FASTERAC_LIBS} ${ROOT_LIBS}

TARGET       = faster2root

# Dossier contenant les fichiers source
SRCDIR = src

# Liste des fichiers source (.cpp) dans le dossier SRCDIR
SOURCES = $(wildcard $(SRCDIR)/*.cpp)

# Liste des fichiers objets correspondants aux fichiers source
OBJECTS = $(SOURCES:.cpp=.o)

# Règle par défaut : génère l'exécutable
all: $(TARGET)

# Règle pour générer l'exécutable
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) -o $(TARGET) $(LIBS)

# Règle générique pour la compilation des fichiers source en fichiers objets
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Nettoyage des fichiers objets et de l'exécutable
clean:
	rm -f $(OBJECTS) $(TARGET)
