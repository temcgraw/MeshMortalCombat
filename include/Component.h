#ifndef COMPONENT_H
#define COMPONENT_H


// forward declaration
class SceneObject;

class IComponent {
public:
    virtual ~IComponent() = default;
    virtual void Tick(float deltaTime) = 0;
    bool active = true;
    SceneObject * getOwner() {
        return owner;
    }
    void setOwner(SceneObject * _owner) {
        owner = _owner;
    }
    bool hasOwner() {
        return owner != nullptr;
    }
protected:
    SceneObject * owner = nullptr;
};

enum renderQueue {
    OPAQUE = 0,
    SKYBOX = 1,
    TRANSPARENT = 2
    
};



class RenderComponent : public IComponent {
public:
    virtual void Render() = 0;
    virtual void Tick(float deltaTime) override {
        return;
    }
    // render priority, higher value means it will be rendered later
    // default value is 0, 1 for skybox, 2 for transparent objects
    renderQueue renderPriority = OPAQUE; 
};

// compute component, used for objects which need to be updated by compute shader
class ComputeComponent : public IComponent {
public:
    virtual void Compute() = 0;
    virtual void Tick(float deltaTime) override {
        return;
    }
};


















#endif

