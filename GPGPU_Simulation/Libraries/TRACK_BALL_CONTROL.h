

#pragma once 
#include <Windows.h>
#include <VECTOR3_T.h>
#include <QUATERNION_T.h>

class TRACK_BALL_CONTROL
{
public:
	
	Vec3T position;
	QuaterT rotation;

	POINT pre_mouse_point;
	POINT cur_mouse_point;

	float dx, dy, dz, dt;

public:

	TRACK_BALL_CONTROL()
	{
		position = Vec3T(0.0f, 0.0f, 0.0f);
		rotation = QuaterT::IDENTITY;

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
			Vec3T angular_vector = Vec3T(1,0,0) * dy * dt + Vec3T(0,1,0) * dx * dt;

			QuaterT q; q = QuaterT::IDENTITY;

			const float magnitude = angular_vector.Magnitude();

			if(magnitude > FLT_EPSILON)
			{
				q.w = cos(magnitude*0.5);

				Vec3T axis = (sinf(0.5*magnitude)/magnitude)*angular_vector;

				q.x = axis.x; q.y = axis.y; q.z = axis.z;
			}

			rotation =  q * rotation; rotation.Normalize();
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