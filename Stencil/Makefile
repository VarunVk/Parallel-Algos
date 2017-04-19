
NVCC=nvcc

CFLAGS += \
		-std=c99 \
		-O3 \
    -funroll-loops \
		-march=native \
    -fopenmp \
		-Wall \
    -Wno-unused-variable

CUDA_FLAGS += \
		-O3 \
    -Wno-deprecated-gpu-targets

LD_FLAGS += -lm -lgomp


# Directories
SRC_DIR = src
OBJ_DIR = obj
EXE_DIR = bin

C_SRC_FILES		= $(wildcard $(SRC_DIR)/*.c)
C_OBJ_FILES = $(addprefix $(OBJ_DIR)/,$(notdir $(C_SRC_FILES:.c=.o)))

CUDA_SRC_FILES = $(wildcard $(SRC_DIR)/*.cu)
CUDA_OBJ_FILES = $(addprefix $(OBJ_DIR)/,$(notdir $(CUDA_SRC_FILES:.cu=.o)))

# Output files
EXE_FILES = stencil

ALL_OBJ_FILES = $(C_OBJ_FILES) $(CUDA_OBJ_FILES)

all: $(EXE_FILES) $(ALL_OBJ_FILES)


dirguard=@mkdir -p $(@D)


$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(dirguard)
	$(CC) $(CFLAGS) -c -o $@ $<

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cu
	$(dirguard)
	$(NVCC) $(CUDA_FLAGS) -c -o $@ $<


# executable creatiion
stencil: $(C_OBJ_FILES) $(CUDA_OBJ_FILES)
	mkdir -p $(EXE_DIR)
	$(NVCC) $(CUDA_FLAGS) $(C_OBJ_FILES) $(CUDA_OBJ_FILES) -o $(EXE_DIR)/stencil $(LD_FLAGS)


clean:
	rm -rf $(OBJ_DIR)

