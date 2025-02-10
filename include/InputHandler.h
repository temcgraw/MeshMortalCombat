// InputHandler.h
#include <functional>
#include <glad/glad.h>
#include <GLFW/glfw3.h>


class InputHandler {
public:
    InputHandler(GLFWwindow* window) {
        this->window = window;
        // let current instance be the user pointer of the window
        // so InputHandler is in charge of the window callbacks
        glfwSetWindowUserPointer(window, this);
        // register callback functions
        glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
        glfwSetCursorPosCallback(window, mouse_callback);
        glfwSetScrollCallback(window, scroll_callback);

        // tell GLFW to capture our mouse
        //glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
        //glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    }

    
    // set callback functions:
    void setMouseMovementCallback(std::function<void(double, double)> callback) {
        this->mouseMovementCallback = callback;
    }
    void setScrollCallback(std::function<void(double)> callback) {
        this->scrollCallback = callback;
    }
    void setFramebufferSizeCallback(std::function<void(int, int)> callback) {
        this->framebufferSizeCallback = callback;
    }
    void setKeyboardActionExecution(std::function<void()> callback) {
        this->keyboardActionExecution = callback;
    }
    void processKeyboardInput() {
        this->keyboardActionExecution();
    }

private:
    GLFWwindow* window;
    std::function<void(float, float)> mouseMovementCallback;
    std::function<void(float)> scrollCallback;
    std::function<void(int, int)> framebufferSizeCallback;
    std::function<void()> keyboardActionExecution;

    static void mouse_callback(GLFWwindow* window, double xpos, double ypos)
    {
	    // get class instance pointer (input handler itself)
	    InputHandler* handler = static_cast<InputHandler*>(glfwGetWindowUserPointer(window));

	    // call its member function (registered from outside, with lambda capture)
        if (handler) {
            handler->mouseMovementCallback(xpos, ypos);
	    }
    }
    static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
    {
        InputHandler* handler = static_cast<InputHandler*>(glfwGetWindowUserPointer(window));
        if (handler) {
            handler->scrollCallback(yoffset);
        }
    }
    static void framebuffer_size_callback(GLFWwindow* window, int width, int height)
    {
        InputHandler* handler = static_cast<InputHandler*>(glfwGetWindowUserPointer(window));
        if (handler) {
            handler->framebufferSizeCallback(width, height);
        }
    }
};
