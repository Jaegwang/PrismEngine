

#pragma once 
#include <Windows.h>
#include "MATH_CORE.h"

class TRACK_BALL_CONTROL
{
public:
	
	TV3 position;
	TQ rotation;

	POINT pre_mouse_point;
	POINT cur_mouse_point;

	float dx, dy, dz, dt;

public:

	TRACK_BALL_CONTROL()
	{
		position = TV3(0.0f, 0.0f, 0.0f);
		rotation = TQ();

		dx = dy = dz = 0.0f;
		dt = 0.005f;
	}
	~TRACK_BALL_CONTROL()
	{}

	void GetInputState()
	{
		if(GetAsyncKeyState(VK_LBUTTON) & 0x8001 || GetAsyncKeyState(VK_RBUTTON) & 0x8001 || GetAsyncKeyState(VK_MBUTTON) & 0x8001)
		{
			pre_mouse_point = cur_mouse_point;
			GetCursorPos(&cur_mouse_point);

			dx = (float)cur_mouse_point.x - (float)pre_mouse_point.x;
			dy = (float)cur_mouse_point.y - (float)pre_mouse_point.y;
		}
		else
		{
			GetCursorPos(&cur_mouse_point); 
			GetCursorPos(&pre_mouse_point);
			dx = dy = 0.0f;
		}

		if(GetAsyncKeyState(VK_MENU) & 0x8001 && GetAsyncKeyState(VK_LBUTTON) & 0x8001)
		{
			TV3 angular_vector = TV3(1,0,0) * dy * dt + TV3(0,1,0) * dx * dt;

			TQ q;

			const float magnitude = glm::length(angular_vector);

			if(magnitude > FLT_EPSILON)
			{
				q.w = cos(magnitude*0.5);

				TV3 axis = (sinf(0.5*magnitude)/magnitude)*angular_vector;

				q.x = axis.x; q.y = axis.y; q.z = axis.z;
			}

			rotation =  q * rotation; 
			rotation = glm::normalize(rotation);
		}

		if(GetAsyncKeyState(VK_MENU) & 0x8001 && GetAsyncKeyState(VK_MBUTTON) & 0x8001)
		{
			position.x += dx * dt;
			position.y -= dy * dt;
		}

		if(GetAsyncKeyState(VK_CONTROL) & 0x8001 && GetAsyncKeyState(VK_LBUTTON) & 0x8001)
		{
			position.z += dy * dt;
		}
	}
};