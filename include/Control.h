#ifndef CONTROL_H
#define CONTROL_H

#include <Camera.h>


// we need a class to have a callback function for mouse movement to update the camera
class CameraController {
public:
    CameraController(Camera* camera)
        : camera(camera), lastX(0.0f), lastY(0.0f), firstMouse(true) {}

    virtual void rotateCamera(double xposIn, double yposIn) {
        float xpos = static_cast<float>(xposIn);
        float ypos = static_cast<float>(yposIn);

        if (firstMouse) {
            lastX = xpos;
            lastY = ypos;
            firstMouse = false;
        }

        float xoffset = xpos - lastX;
        float yoffset = lastY - ypos;  // reversed since y-coordinates go from bottom to top

        lastX = xpos;
        lastY = ypos;
        if (enableCameraControl) {
            camera->ProcessMouseMovement(xoffset, yoffset);
        }
    }

    virtual void adjustfov(float yoffset) {
        if (enableCameraControl){
            camera->ProcessMouseScroll(static_cast<float>(yoffset));
        }
    }

    // just warp the camera's processKeyboard function for better hierarchy and readability
    virtual void moveCamera(float deltaTime, enum Camera_Movement direction) {
        if (enableCameraControl) {
            camera->ProcessKeyboard(direction, deltaTime);
        }
    }


    void setEnableCameraControl(bool enable) {
        enableCameraControl = enable;
    }

    bool canControlCamera() {
        return enableCameraControl;
    }

protected:
    Camera* camera;
    float lastX, lastY;
    bool firstMouse;
    bool enableCameraControl = true;
};



void keyboardActions(ExampleStateMachine * state_machine,CameraController * camera_controller , float deltaTime, GLFWwindow* window
)
{    
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS){
        glfwSetWindowShouldClose(window, true);
    }
    std::string current_state = state_machine->get_current_state()->name;
    if (current_state != "CPU ray-tracing") { // while in CPU ray tracing state, we disable the keyboard input
        if (camera_controller->canControlCamera()) {
            bool camera_moved = false;
            if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS){
                camera_controller->moveCamera(deltaTime, FORWARD);
                camera_moved = true;
            }

            if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS){
                camera_controller->moveCamera(deltaTime, BACKWARD);
                camera_moved = true;
            }
            if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS){
                camera_controller->moveCamera(deltaTime, LEFT);
                camera_moved = true;
            }
            if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS){
                camera_controller->moveCamera(deltaTime, RIGHT);
                camera_moved = true;
            }
            if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS){
                camera_controller->moveCamera(deltaTime, UP);
                camera_moved = true;
            }
            if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS){
                camera_controller->moveCamera(deltaTime, DOWN);
                camera_moved = true;
            }
        }

        if (glfwGetKey(window, GLFW_KEY_F) == GLFW_PRESS) {
            camera_controller->setEnableCameraControl(false);
            glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
	    }
        if (glfwGetKey(window, GLFW_KEY_G) == GLFW_PRESS) {
            camera_controller->setEnableCameraControl(true);
            glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
        }

        if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS) {
            state_machine->request_state2();
        }
        if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS) {
            state_machine->request_state1();
        }
        if (glfwGetKey(window, GLFW_KEY_3) == GLFW_PRESS) {
            state_machine->request_state3();
        }
    }
}







#endif
