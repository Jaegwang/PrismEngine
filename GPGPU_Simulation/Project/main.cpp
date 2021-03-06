
#include "READ_WRITE_PNG.h"

#include <GL\glut.h>
#include <iostream>
#include <amp.h>

#include "PARTICLE_MANAGER_3D.h"
#include "PARTICLE.h"
#include "MPM_FLUID_SOLVER.h"
#include "CAPTURE_MANAGER.h"

#include "MATH_CORE.h"
#include "GRID.h"
#include "TRACK_BALL_CONTROL.h"

#include "GRID_DATA.h"
#include "FIELD_ENCODED.h"
#include "PROJECTION_METHOD.h"
#include "STABLE_FLUID_SOLVER.h"

#include "SPARSE_MATRIX.h"


using namespace concurrency;


TRACK_BALL_CONTROL track_ball;
MPM_FLUID_SOLVER mpm_solver;

CAPTURE_MANAGER capture_manager;

FIELD_ENCODED<TS> field;

STABLE_FLUID_SOLVER stable_fluid;


static const int window_w = 800;
static const int window_h = 600;


static bool is_playing      = false;
static bool is_capture      = true;
static bool is_capture_flag = false;


void display();
void idle();
void reshape(int w, int h);
void mouseButton(int button,int state,int x,int y);
void mouseMotion(int x,int y);
void keyboard(unsigned char key, int x, int y);
void init();
void light();


int main(int argc, char **argv)
{
//	omp_set_dynamic(0);
//	omp_set_num_threads(4);
//	int num = omp_get_max_threads();

	TV3 min0(0,0,0);
	TV3 max0(1,1,1);

	std::string path = "no";

	GRID grid;
	grid.Initialize(min0, max0, 100, 100, 100, 2);

	stable_fluid.Initialize(grid);

	stable_fluid.SeedParticlesFromSphere(TV3(0.5, 0.9, 0.5), 0.1, TV3(0.0, -2, 0.0), 20000);
//	stable_fluid.SeedParticlesFromSphere(TV3(0.7, 0.5, 0.7), 0.1, TV3(0.0, -2, 0.0), 20000);

	stable_fluid.SeedParticlesFromHeight(0.3, 3);

	mpm_solver.Initialize(min0, max0, 100, 100, 100, 2, 5000000);
	capture_manager.Initialize(path);

	
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(window_w, window_h);
	glutCreateWindow("Stable Fluids");

	glutReshapeFunc(reshape);
	glutIdleFunc(idle);
	glutDisplayFunc(display);

	glutMouseFunc(mouseButton);
	glutMotionFunc(mouseMotion);
	glutKeyboardFunc(keyboard);

	init();
	light();

	glutMainLoop();
}

void init()
{

}

void display()
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	const TV3& world_center = (mpm_solver.grid_.max_ + mpm_solver.grid_.min_)*(TS)0.5;

	glTranslatef(0, 0, -2.0f);
	
	glTranslatef(track_ball.position.x, track_ball.position.y, track_ball.position.z);	

	glRotatef((float)acos(track_ball.rotation.w)*(float)360/(float)3.141592, track_ball.rotation.x, track_ball.rotation.y, track_ball.rotation.z);	
	glTranslatef(-world_center.x, -world_center.y, -world_center.z);


	glClearColor(0.3, 0.3, 0.3, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

//	glutSolidTeapot(0.5);
//	world_grid.RenderGrid();
//	world_grid.RenderCells();
//	particle_manager.Rendering();

//	mpm_solver.particle_manager_.Rendering();
//	mpm_solver.grid_.RenderGrid();
//	mpm_solver.particle_world_.Render();

//	stable_fluid.Render();

	stable_fluid.Grid().RenderGrid();
	stable_fluid.RenderParticles();
//	stable_fluid.RenderVelocity();
	

	FIELD_ENCODED<TS>* df = (FIELD_ENCODED<TS>*)stable_fluid.density_field_;

//	df->Render();

//	field.Render();

	// capture image and video
	if(is_capture && is_capture_flag)
	{
		static int __capture_frame = 0;
		GLint viewport[4];
		glGetIntegerv(GL_VIEWPORT, viewport);
		capture_manager.CaptureImage(__capture_frame++, viewport[2] - viewport[0], viewport[3] - viewport[1]);
		is_capture_flag = false;
	}

	glFlush();
	glutSwapBuffers();
}

void idle()
{
	track_ball.GetInputState();

	if (is_playing == true)
	{
//		mpm_solver.AdvanceTimeStep((TS)0.01, 1);

//		stable_fluid.SourceDensityFromSphere(TV3(0.5, 0.2, 0.5), 0.05, 1.0, TV3(0.0, 2, 0.0));
	
//		stable_fluid.SeedParticlesFromSphere(TV3(0.5, 0.8, 0.5), 0.1, TV3(0.0, -2, 0.0), 5000);
		stable_fluid.AdvanceOneTimeStepFLIP(0.01);
	//	stable_fluid.AdvanceOneTimeStep(0.01);

		is_capture_flag = true;
	}

	glutPostRedisplay();
}

void reshape(int w, int h)
{
	glViewport(0, 0, w, h);

	// projection
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(60.0f, (double) w / h, 0.1, 1000);
}

void mouseButton(int button,int state,int x,int y)
{}

void mouseMotion(int x,int y)
{}

void keyboard(unsigned char key, int x, int y)
{
	key = tolower(key);

	switch (key)
	{
		case 'q': exit(0);	                    break;
		case 'p': is_playing = !is_playing;		break;
		case 'c': is_capture = !is_capture;		break;
		case 'v': capture_manager.MakeVideo();  break;
	}

	glutPostRedisplay();
}

void light()
{
	glEnable(GL_SMOOTH);
	glEnable(GL_NORMALIZE);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);

	// default light

	GLfloat gl_light_position[]={0,0,0,1};
	GLfloat gl_light_ambient[] ={1,1,1,1};
	GLfloat gl_light_diffuse[] ={1,1,1,1};
	GLfloat gl_light_specular[]={0,0,0,0};

	glLightfv(GL_LIGHT0,GL_POSITION,gl_light_position);
	glLightfv(GL_LIGHT0,GL_AMBIENT,gl_light_ambient);			
	glLightfv(GL_LIGHT0,GL_DIFFUSE,gl_light_diffuse);			
	glLightfv(GL_LIGHT0,GL_SPECULAR,gl_light_specular);			
	glEnable(GL_LIGHT0);

	// default material

	GLfloat gl_material_ambient[]={0.2,0.2,0.2,1};
	GLfloat gl_material_diffuse[]={0.3,0.3,0.3,1}; 
	GLfloat gl_material_specular[]={1,1,1};
	GLfloat gl_material_shininess=100;
	GLfloat gl_material_emission[]={0,0,0};

	glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);			//Set Material properties
	glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,gl_material_ambient);		//Set Material Ambient
	glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,gl_material_specular);	//Set Material Specular
	glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,gl_material_diffuse);		//Set Material Diffuse
	glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,gl_material_shininess);	//Set Material Shininess
	glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,gl_material_emission);	//Set Material Emission
}