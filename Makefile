# **************************************************************************** #
#                                                                              #
#                                                                              #
#    Makefile                                          Personal Website        #
#                                                     ##################       #
#    By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||       #
#                                                     ##################       #
#    Created: 2022/11/09 19:28:22 by Zian Huang                                #
#    Updated: 2022/11/16 16:22:42 by Zian Huang                                #
#                                                                              #
# **************************************************************************** #

src_dir = src
obj_dir = obj
bin_dir = bin
compiler = clang++
compile_std = -std=c++11

objects = $(obj_dir)/flux_func.o $(obj_dir)/slic_2d_euler_solver.o $(obj_dir)/vec_transform.o $(obj_dir)/exec_interface.o

all: flux_func.o slic_2d_euler_solver.o vec_transform.o exec_interface.o
	$(compiler) $(compile_std) $(objects) -o $(bin_dir)/slic_solver_2d

exec_interface.o: $(src_dir)/exec_interface.cc
	$(compiler) -c $(compile_std) -O2 $(src_dir)/exec_interface.cc -o $(obj_dir)/exec_interface.o

flux_func.o: $(src_dir)/flux_func.cc
	$(compiler) -c $(compile_std) -O2 $(src_dir)/flux_func.cc -o $(obj_dir)/flux_func.o

slic_2d_euler_solver.o: $(src_dir)/slic_2d_euler_solver.cc
	$(compiler) -c $(compile_std) -O2 $(src_dir)/slic_2d_euler_solver.cc -o $(obj_dir)/slic_2d_euler_solver.o

vec_transform.o: $(src_dir)/vec_transform.cc
	$(compiler) -c $(compile_std) -O2 $(src_dir)/vec_transform.cc -o $(obj_dir)/vec_transform.o

clean:
	rm $(bin_dir)/slic_solver_2d $(objects)
