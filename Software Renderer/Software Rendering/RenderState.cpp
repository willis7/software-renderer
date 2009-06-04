#include "RenderState.h"

RenderState::RenderState(void)
{
	eyePos = Vector3(0.0, 0.0, 0.0);
	eyeDir = Vector3(0.0, 0.0, 1.0);
	depthTest = false;
	cullMode  = CullNone;
	lightingEnabled  = false;
	texturingEnabled = false;
	style = Flat;
}


