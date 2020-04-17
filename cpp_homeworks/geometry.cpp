#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
 
using namespace std;

struct Point {
    double x;
    double y;
    Point() : x(0), y(0) {};
    Point(const double& x, const double& y) : x(x), y(y) {};
    bool operator==(const Point& other) const;
    bool operator!=(const Point& other) const;
};

bool Point::operator==(const Point& other) const {
    return x == other.x && y == other.y;
}

bool Point::operator!=(const Point& other) const {
    return !(*this == other);
}

double distance(const Point& left, const Point& right) {
    return sqrt((left.x - right.x) * (left.x - right.x) + (left.y - right.y) * (left.y - right.y));
}

void rotatePoint(Point& point, const Point& center, double& angle) {
    point.x -= center.x;
    point.y -= center.y;
    
    point.x = cos(angle) * point.x - sin(angle) * point.y;
    point.y = sin(angle) * point.x + cos(angle) * point.y;

    point.x += center.x;
    point.y += center.y;
}

void scalePoint(Point& point, const Point& center, double coef) {
    Point homo_vector(point.x - center.x, point.y - center.y);
    homo_vector.x *= coef;
    homo_vector.y *= coef;
    point = Point(center.x + homo_vector.x, center.y + homo_vector.y);
}

class Line {
public:
    Line() : a(0), b(0), c(0) {};
    Line(const Point& first, const Point& second);
    Line(const double& angle_coef, const double& shift);
    Line(const Point& first, const double& angle_coef);
    double get_coef() const;
    bool operator==(const Line& other) const;
    friend Point lineIntersection(const Line& first, const Line& second);

private:
    // Ax + By + C = 0
    double a, b, c;
};

Line::Line(const Point& first, const Point& second) {
    a = 1.0;
    b = (second.x - first.x) / (second.y - first.y);
    c = (first.x * second.y - first.y * second.x) / (second.y - first.y);
}

Line::Line(const double& angle_coef, const double& shift) {
    Point first(shift, 0);
    Point second(shift + 1, angle_coef);
    a = 1.0;
    b = (second.x - first.x) / (second.y - first.y);
    c = (first.x * second.y - first.y * second.x) / (second.y - first.y);
}

Line::Line(const Point& first, const double& angle_coef) {
    Point second(first.x + 1, first.y + angle_coef);
    a = 1.0;
    b = (second.x - first.x) / (second.y - first.y);
    c = (first.x * second.y - first.y * second.x) / (second.y - first.y);
}

double Line::get_coef() const {
    return -(a / b);
}

bool Line::operator==(const Line& other) const {
    return a == other.a && b == other.b && c == other.c;
}

Point lineIntersection(const Line& first, const Line& second) {
    double x = (first.b * second.c - first.c * second.b) / (first.a * second.b - first.b * second.a);
    double y = (-first.a * x - first.c) / first.b;
    return Point(x, y);
}

void reflexPoint(Point& point, const Line& line) {
    Line ortholine = Line(point, 1 / line.get_coef());
    Point intersection = lineIntersection(line, ortholine);
    double add_x = intersection.x - point.x;
    double add_y = intersection.y - point.y;
    point.x += 2 * add_x;
    point.y += 2 * add_y;
}

void reflexPoint(Point& point, const Point& center) {
    double add_x = center.x - point.x;
    double add_y = center.y - point.y;
    point.x += 2 * add_x;
    point.y += 2 * add_y;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum ShapeTypes {
    ST_POLYGON,
    ST_ELLIPSE
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Shape {
public:
    virtual double perimeter() const = 0;
    virtual double area() const = 0;
    virtual bool operator==(const Shape& other) const = 0;
    virtual bool isCongruentTo(const Shape& another) const = 0;
    virtual bool isSimilarTo(const Shape& another) const = 0;
    virtual bool containsPoint(Point point) const = 0;

    virtual void rotate(const Point& center, double& angle) = 0;
    virtual void reflex(const Point& center) = 0;
    virtual void reflex(Line axis) = 0;
    virtual void scale(Point center, double coefficient) = 0;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Ellipse : public Shape {
public:
    Ellipse() {};
    Ellipse(Point foc1, Point foc2, double _dist);
    pair<Point, Point> focuses() const;
    pair<Line, Line> directrices() const;
    double eccentricity() const;
    Point center() const;

    double perimeter() const override;
    double area() const override;
    bool operator==(const Shape& other) const override;
    bool isCongruentTo(const Shape& another) const override;
    bool isSimilarTo(const Shape& another) const override;
    bool containsPoint(Point point) const override;

    void rotate(const Point& center, double& angle) override;
    void reflex(const Point& center) override;
    void reflex(Line axis) override;
    void scale(Point center, double coefficient) override;

private:
    pair<Point, Point> focus;
    double dist;
    double a() const;
    double b() const;
    double c() const;
};

Ellipse::Ellipse(Point foc1, Point foc2, double _dist) {
    focus = make_pair(foc1, foc2);
    dist = _dist;
}

pair<Point, Point> Ellipse::focuses() const {
    return focus;
}

double Ellipse::a() const {
    return dist / 2;
}

double Ellipse::b() const {
    return a() * a() - c() * c();
}

double Ellipse::c() const {
    return distance(focus.first, focus.second);
}

double Ellipse::eccentricity() const {
    return c() / a();
}

double Ellipse::perimeter() const {
    return (M_PI * a() * b() + a() - b()) / (a() + b());
}

double Ellipse::area() const {
    return M_PI * a() * b();
}

void Ellipse::rotate(const Point& center, double& angle) {
    rotatePoint(focus.first, center, angle);
    rotatePoint(focus.second, center, angle);
}

void Ellipse::reflex(const Point& center) {
    reflexPoint(focus.first, center);
    reflexPoint(focus.second, center);
}

void Ellipse::reflex(Line axis) {
    reflexPoint(focus.first, axis);
    reflexPoint(focus.second, axis);
}

void Ellipse::scale(Point center, double coefficient) {
    scalePoint(focus.first, center, coefficient);
    scalePoint(focus.second, center, coefficient);
    dist *= coefficient;
}

///////////////////////////////////////////////////////////////////
                                                                ///
pair<Line, Line> Ellipse::directrices() const {                 ///
    double coef = Line(focus.first, focus.second).get_coef();   ///
    double e = eccentricity();                                  ///
    coef = -(1 / coef);                                         ///
    //Point p1(center().x + );                                  ///
}                                                               ///
///////////////////////////////////////////////////////////////////

Point Ellipse::center() const {
    return Point((focus.first.x + focus.second.x) / 2, (focus.first.y + focus.second.y) / 2);
}

class Circle : public Ellipse {
public:
    Circle(const Point& center, double rad) : center(center), rad(rad), Ellipse(center, center, 2 * rad) {};
    double radius() const;

    double perimeter() const override;
    double area() const override;

    void rotate(const Point& _center, double& angle) override;
    void reflex(const Point& _center) override;
    void reflex(Line axis) override;
    void scale(Point _center, double coefficient) override;

private:
    Point center;
    double rad;
};

void Circle::rotate(const Point& _center, double& angle) {
    Point new_center = center;
    rotatePoint(new_center, _center, angle);
    *this = Circle(new_center, rad);
}

void Circle::reflex(const Point& _center) {
    Point new_center = center;
    reflexPoint(new_center, _center);
    *this = Circle(new_center, rad);
}

void Circle::reflex(Line axis) {
    Point new_center = center;
    reflexPoint(new_center, axis);
    *this = Circle(new_center, rad);
}

void Circle::scale(Point _center, double coefficient) {
    Point new_center = center;
    scalePoint(new_center, _center, coefficient);
    *this = Circle(new_center, rad);
}

double Circle::radius() const {
    return rad;
}

double Circle::perimeter() const {
    return 2 * M_PI * rad;
}

double Circle::area() const {
    return M_PI * rad * rad;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Polygon : public Shape {
public:
    Polygon(const vector<Point>& points);
    template<typename... Args>
    Polygon(Args...);
    size_t verticesCount() const;
    vector<Point> getVertices() const;
    bool isConvex() const;

    double perimeter() const override;
    double area() const override;
    bool operator==(const Shape& other) const override;
    bool isCongruentTo(const Shape& another) const override;
    bool isSimilarTo(const Shape& another) const override;
    bool containsPoint(Point point) const override;

    void rotate(const Point& _center, double& angle) override;
    void reflex(const Point& _center) override;
    void reflex(Line axis) override;
    void scale(Point _center, double coefficient) override;

private:
    vector<Point> vertices;
    void add_vertexes() {};
    template<typename Head, typename... Tail>
    void add_vertexes(const Head& head, const Tail&... tail);
};

Polygon::Polygon(const vector<Point>& points) {
    vertices = points;    
}

template<typename... Args>
Polygon::Polygon(Args... args) {
    add_vertexes(args...);
}

size_t Polygon::verticesCount() const {
    return vertices.size();
}

vector<Point> Polygon::getVertices() const {
    return vertices;
}

bool Polygon::isConvex() const {
    bool flag_neg = false;
    bool flag_pos = false;
    for (size_t i = 0; i < verticesCount(); ++i) {
        size_t j = (i + 1) % verticesCount();
        size_t k = (i + 2) % verticesCount();
        double z = (vertices[j].x - vertices[i].x) * (vertices[j].y - vertices[i].y) 
                - (vertices[k].x - vertices[j].x) * (vertices[k].y - vertices[j].y);
        if (z < 0)
            flag_neg |= true;
        else if (z > 0)
            flag_pos |= true;
        if (flag_neg && flag_pos)
            return false;
    }
    if (flag_neg || flag_pos)
        return true;
    else 
        return false;
}

void Polygon::rotate(const Point& _center, double& angle) {
    for (size_t i = 0; i < vertices.size(); ++i)
        rotatePoint(vertices[i], _center, angle);
}

void Polygon::reflex(const Point& _center) {
    for (size_t i = 0; i < vertices.size(); ++i)
        reflexPoint(vertices[i], _center);
}

void Polygon::reflex(Line axis) {
    for (size_t i = 0; i < vertices.size(); ++i)
        reflexPoint(vertices[i], axis);
}

void Polygon::scale(Point _center, double coef) {
    for (size_t i = 0; i < vertices.size(); ++i)
        scalePoint(vertices[i], _center, coef);
}

template<typename Head, typename... Tail>
void Polygon::add_vertexes(const Head& head, const Tail&... tail) {
    vertices.push_back(head);
    add_vertexes(tail...);
}

double Polygon::perimeter() const {
    double res = distance(vertices[0], vertices[vertices.size() - 1]);
    for (size_t i = 0; i < vertices.size() - 1; ++i)
        res += distance(vertices[i], vertices[i + 1]);
    return res;
}

double Polygon::area() const {
    double res = (vertices[0].x + vertices[vertices.size() - 1].x) * (vertices[0].y - vertices[vertices.size() - 1].y);
    for (size_t i = 0; i < vertices.size() - 1; ++i)
        res += (vertices[i].x + vertices[i + 1].x) * (vertices[i].y - vertices[i + 1].y);
    res /= 2;
    return res;
}

class Triangle : public Polygon {
public:
    Triangle(Point a, Point b, Point c) : Polygon(a, b, c) {};
    Circle circumscribedCircle() const;
    Circle inscribedCircle() const;

    Point centroid() const;
    Point orthocenter() const;

    Line EulerLine() const;
    Circle ninePointsCircle() const;
    
    double perimeter() const override;
    double area() const override;
    bool operator==(const Shape& other) const override;
    bool isCongruentTo(const Shape& another) const override;
    bool isSimilarTo(const Shape& another) const override;
    bool containsPoint(Point point) const override;

    void rotate(const Point& center, double& angle) override;
    void reflex(const Point& center) override;
    void reflex(Line axis) override;
    void scale(Point center, double coefficient) override;
};

Point Triangle::centroid() const {
    vector<Point> verts = getVertices();
    Point mid1 = Point((verts[0].x + verts[1].x) / 2, (verts[0].y + verts[1].y) / 2);
    Point mid2 = Point((verts[1].x + verts[2].x) / 2, (verts[1].y + verts[2].y) / 2);
    Line median1 = Line(verts[2], mid1);
    Line median2 = Line(verts[0], mid2);
    return lineIntersection(median1, median2);
}

Point Triangle::orthocenter() const {
    vector<Point> verts = getVertices();
    Line edge1 = Line(verts[1], verts[2]);
    Line edge2 = Line(verts[0], verts[2]);
    Line height1(verts[0], -edge1.get_coef());
    Line height2(verts[1], -edge2.get_coef());
    return lineIntersection(height1, height2);
}

class Rectangle : public Polygon {
public:
    Rectangle(Point x, Point y, double length, double width);
    Rectangle(Point a, Point b, Point c, Point d) : Polygon(a, b, c, d) {};
    Point center() const;
    pair<Line, Line> diagonals() const;

    double perimeter() const override;
    double area() const override;
    bool operator==(const Shape& other) const override;
    bool isCongruentTo(const Shape& another) const override;
    bool isSimilarTo(const Shape& another) const override;
    bool containsPoint(Point point) const override;

    void rotate(const Point& center, double& angle) override;
    void reflex(const Point& center) override;
    void reflex(Line axis) override;
    void scale(Point center, double coefficient) override;

private:

};

Rectangle::Rectangle(Point x, Point y, double length, double width) {
    if (length > width)
        swap(length, width);
    double k = distance(x, y) / sqrt(width * width + length * length);
    length *= k;
    width *= k;
    Point center = x.y < y.y ? x : y;
    Line diagonal(x, y);
    Line unrotated_diag;
    Point a, b;
    if (x.y < y.y) {
        a = Point(x.x + length, x.y);
        b = Point(x.x, x.y + width);
        unrotated_diag = Line(x, Point(x.x + length, x.y + width));
        double angle = (diagonal.get_coef - unrotated_diag.get_coef) / (1 + diagonal.get_coef * unrotated_diag.get_coef);
        rotatePoint(a, x, angle);
        rotatePoint(b, x, angle);
    } else {
        a = Point(y.x + width, y.y);
        b = Point(y.x, y.y + length);
        unrotated_diag = Line(x, Point(y.x + width, y.y + width));
        double angle = (diagonal.get_coef - unrotated_diag.get_coef) / (1 + diagonal.get_coef * unrotated_diag.get_coef);
        rotatePoint(a, y, angle);
        rotatePoint(b, y, angle);
    }
    *this = Rectangle(x, a, y, b);
}

pair<Line, Line> Rectangle::diagonals() const {
    return make_pair(Line(getVertices()[0], getVertices()[2]), Line(getVertices()[1], getVertices()[3]));
}

Point Rectangle::center() const {
    return lineIntersection(diagonals().first, diagonals().second);
}

class Square : public Rectangle {
public:
    Square(Point a, Point b) : Rectangle(a, b, 1, 1) {};
    Circle circumscribedCircle() const;
    Circle inscribedCircle() const;
};

Circle Square::inscribedCircle() const {
    return Circle(center(), distance(getVertices()[0], getVertices()[1]) / 2);
}

Circle Square::circumscribedCircle() const {
    return Circle(center(), distance(center(), getVertices()[0]));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main() {
    Polygon a(Point(), Point(), Point(), Point());
    system("pause");
    return 0;
}