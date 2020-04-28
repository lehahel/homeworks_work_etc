#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
 
using namespace std;

#define EPSILON 1e-10;

double toRad(double angle) {
    return (angle / 180) * M_PI;
}

bool areSame(const double& first, const double& second) {
    return fabs(first - second) < EPSILON;
}

struct Point {
    double x;
    double y;
    Point() : x(0), y(0) {};
    Point(const double& x, const double& y) : x(x), y(y) {};
    bool operator==(const Point& other) const;
    bool operator!=(const Point& other) const;
};

bool Point::operator==(const Point& other) const {
    return areSame(x, other.x) && areSame(y, other.y);
}

bool Point::operator!=(const Point& other) const {
    return !(*this == other);
}

double bandistance(const Point& left, const Point& right) {
    return sqrt((left.x - right.x) * (left.x - right.x) + (left.y - right.y) * (left.y - right.y));
}

void rotatePoint(Point& point, const Point& center, double& angle) {
    point.x -= center.x;
    point.y -= center.y;
    
    double x = cos(toRad(angle)) * point.x - sin(toRad(angle)) * point.y;
    double y = sin(toRad(angle)) * point.x + cos(toRad(angle)) * point.y;

    point.x = x;
    point.y = y;

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
    Line(double angle_coef, double shift);
    Line(const Point& first, double angle_coef);
    double get_coef() const;
    bool operator==(const Line& other) const;
    friend Point lineIntersection(const Line& first, const Line& second);
    double a, b, c;
};

Line::Line(const Point& first, const Point& second) {
    a = first.y - second.y;
    b = second.x - first.x;
    c = first.x * second.y - second.x * first.y;
}

Line::Line(double angle_coef, double shift) {
    Point first(0, shift);
    Point second(1, shift + angle_coef);
    *this = Line(first, second);
}

Line::Line(const Point& first, double angle_coef) {
    Point second(first.x + 1, first.y + angle_coef);
    *this = Line(first, second);
}

double Line::get_coef() const {
    return -(a / b);
}

bool Line::operator==(const Line& other) const {
    if (areSame(a, 0) && !areSame(other.a, 0))
        return false;
    if (!areSame(a, 0) && areSame(other.a, 0))
        return false;
    if (areSame(b, 0) && !areSame(other.b, 0))
        return false;
    if (!areSame(b, 0) && areSame(other.b, 0))
        return false;
    if (areSame(c, 0) && !areSame(other.c, 0))
        return false;
    if (!areSame(c, 0) && areSame(other.a, 0))
        return false;
    if (areSame(a, 0))
        return areSame(-c / b, -other.c / other.b);
    if (areSame(b, 0))
        return areSame(-c / a, -other.c / other.a);
    if (areSame(c, 0))
        return areSame(a / other.a, b / other.b);
    return areSame(a / other.a, b / other.b) && areSame(b / other.b, c / other.c);
}

Point lineIntersection(const Line& first, const Line& second) {
    double x = (first.b * second.c - first.c * second.b) / (first.a * second.b - first.b * second.a);
    double y = (-first.a * x - first.c) / first.b;
    return Point(x, y);
}

void reflexPoint(Point& point, const Line& line) {
    Line ortholine = Line(point, -1 / line.get_coef());
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
    bool operator!=(const Shape& other) const;

    virtual void rotate(const Point& center, double angle) = 0;
    virtual void reflex(const Point& center) = 0;
    virtual void reflex(Line axis) = 0;
    virtual void scale(Point center, double coefficient) = 0;

    ShapeTypes type;
    virtual ~Shape() {};
};

bool Shape::operator!=(const Shape& other) const {
    return !(*this == other);
}

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

    void rotate(const Point& center, double angle) override;
    void reflex(const Point& center) override;
    void reflex(Line axis) override;
    void scale(Point center, double coefficient) override;

    double focal_distance() const;

    ~Ellipse() override {};

private:
    pair<Point, Point> focus;
    double dist;
    double a() const;
    double b() const;
    double c() const;
};

double Ellipse::focal_distance() const {
    return dist;
}

Ellipse::Ellipse(Point foc1, Point foc2, double _dist) {
    focus = make_pair(foc1, foc2);
    dist = _dist;
    type = ST_ELLIPSE;
}

bool Ellipse::isSimilarTo(const Shape& other) const {
    if (other.type != ST_ELLIPSE)
        return false;
    Ellipse shape = dynamic_cast<const Ellipse&>(other);

    return a() / b() == shape.a() / shape.b() || b() / a() == shape.a() / shape.b();
}

bool Ellipse::isCongruentTo(const Shape& another) const {
    return isSimilarTo(another) && area() == another.area();
}

bool Ellipse::operator==(const Shape& other) const {
    if (other.type != ST_ELLIPSE)
        return false;
    Ellipse shape = dynamic_cast<const Ellipse&>(other);
    if ((focus.first == shape.focus.first && focus.second == shape.focus.second 
    || focus.first == shape.focus.second && focus.second == shape.focus.first) && areSame(dist, shape.dist))
        return true;
    return false;
}

pair<Point, Point> Ellipse::focuses() const {
    return focus;
}

double Ellipse::a() const {
    return dist / 2;
}

double Ellipse::b() const {
    return sqrt(a() * a() - c() * c());
}

double Ellipse::c() const {
    return bandistance(focus.first, focus.second) / 2;
}

double Ellipse::eccentricity() const {
    return c() / a();
}

double Ellipse::perimeter() const {
    return M_PI * (3 * (a() + b()) - sqrt(3. * a() * a() + 10. * a() * b() + 3. * b() * b()));
}

double Ellipse::area() const {
    return M_PI * a() * b();
}

void Ellipse::rotate(const Point& center, double angle) {
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
                                                                
pair<Line, Line> Ellipse::directrices() const {
    double e = eccentricity();
    Point cent = center();
    if (focus.first.x == focus.second.x) {
        Point first(cent.x, cent.y - a() / e);
        Point second(cent.x, cent.y + a() / e);
        return make_pair(Line(first, 0), Line(second, 0));
    }
    if (focus.first.y == focus.second.y) {
        Point first_low(cent.x - a() / e, cent.y);
        Point first_high(cent.x - a() / e, cent.y + 1);
        Point second_low(cent.x + a() /e, cent.y);
        Point second_high(cent.x + a() / e, cent.y + 1);
        return make_pair(Line(first_low, first_high), Line(second_low, second_high));
    }
    double coef = Line(focus.first, focus.second).get_coef();
    double x = 1;
    double y = coef;
    double k = ((a() / e) * (a() / e)) / (x * x + y * y);
    x *= k;
    y *= k;
    Point dir_l(center().x - x, center().y - y);
    Point dir_r(center().x + x, center().y + y);
    return make_pair(Line(dir_l, -1 / coef), Line(dir_r, -1 / coef));
}                                                               

Point Ellipse::center() const {
    return Point((focus.first.x + focus.second.x) / 2, (focus.first.y + focus.second.y) / 2);
}

bool Ellipse::containsPoint(Point point) const {
    return bandistance(point, focus.first) + bandistance(point, focus.second) < dist;
};

class Circle : public Ellipse {
public:
    Circle(const Point& center, double _rad) : Ellipse(center, center, 2 * _rad) {};
    double radius() const;

    double perimeter() const override;
    double area() const override;

    void rotate(const Point& _center, double angle) override;
    void reflex(const Point& _center) override;
    void reflex(Line axis) override;
    void scale(Point _center, double coefficient) override;
};

void Circle::rotate(const Point& _center, double angle) {
    Point new_center = focuses().first;
    rotatePoint(new_center, _center, angle);
    *this = Circle(new_center, focal_distance() / 2);
}

void Circle::reflex(const Point& _center) {
    Point new_center = focuses().first;
    reflexPoint(new_center, _center);
    *this = Circle(new_center, focal_distance() / 2);
}

void Circle::reflex(Line axis) {
    Point new_center = focuses().first;
    reflexPoint(new_center, axis);
    *this = Circle(new_center, focal_distance() / 2);
}

void Circle::scale(Point _center, double coefficient) {
    Point new_center = focuses().first;
    scalePoint(new_center, _center, coefficient);
    *this = Circle(new_center, focal_distance() / 2);
}

double Circle::radius() const {
    return focal_distance() / 2;
}

double Circle::perimeter() const {
    return M_PI * focal_distance();
}

double Circle::area() const {
    return M_PI * (focal_distance() / 2) * (focal_distance() / 2);
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

    void rotate(const Point& _center, double angle) override;
    void reflex(const Point& _center) override;
    void reflex(Line axis) override;
    void scale(Point _center, double coefficient) override;

    ~Polygon() override {};

private:
    vector<Point> vertices;
    void add_vertexes();
    template<typename Head, typename... Tail>
    void add_vertexes(const Head& head, const Tail&... tail);
};

Polygon::Polygon(const vector<Point>& points) {
    vertices = points;
    type = ST_POLYGON;
}

template<typename... Args>
Polygon::Polygon(Args... args) {
    add_vertexes(args...);
    type = ST_POLYGON;
}

size_t Polygon::verticesCount() const {
    return vertices.size();
}

vector<Point> Polygon::getVertices() const {
    return vertices;
}

bool Polygon::isSimilarTo(const Shape& other) const {
    if (other.type != ST_POLYGON)
        return false;
    const Polygon shape = dynamic_cast<const Polygon&>(other);
    vector<Point> vert = getVertices();
    vector<Point> vert_other = shape.getVertices();
    if (vert.size() != vert_other.size())
        return false;
    vector<double> edges;
    vector<double> edges_other;
    for (size_t i = 0; i < vert.size(); ++i) {
        edges.push_back(bandistance(vert[i], vert[(i + 1) % vert.size()]));
        edges_other.push_back(bandistance(vert_other[i], vert_other[(i + 1) % vert.size()]));
    }
    for (size_t i = 0; i < edges.size(); ++i) {
        double k = edges[i] / edges_other[0];
        bool res = true;
        for (size_t j = 0; j < edges.size(); ++j) {
            if (!areSame(edges[(j + i) % edges.size()] / edges_other[j], k)) {
                res = false;
                break;
            }
        }
        if (res)
            return true;
    }
    for (size_t i = 0; i < vert.size(); ++i) {
        double k = edges[i] / edges_other[0];
        bool res = true;
        for (size_t j = 0; j < edges.size(); ++j) {
            if (!areSame(edges[(edges.size() + i - j) % edges.size()] / edges_other[j], k)) {
                res = false;
                break;
            }
        }
        if (res)
            return true;
    }
    return false;
}

bool Polygon::isCongruentTo(const Shape& another) const {
    return isSimilarTo(another) && area() == another.area();
}

bool Polygon::operator==(const Shape& other) const {
    if (other.type != ST_POLYGON)
        return false;
    const Polygon shape = dynamic_cast<const Polygon&>(other);
    vector<Point> vert = getVertices();
    vector<Point> vert_other = shape.getVertices();
    if (vert.size() != vert_other.size())
        return false;
    size_t start_vertex = vert.size();
    for (size_t i = 0; i < vert_other.size(); ++i)
        if (vert[i] == vert_other[0]) {
            start_vertex = i;
            break;
        }
    if (start_vertex == vert.size())
        return false;
    if (vert[(start_vertex + 1) % vert.size()] == vert_other[1]) {
        for (size_t i = 0; i < vert.size(); ++i)
            if (vert[(start_vertex + i) % vert.size()] != vert_other[i])
                return false;
        return true;
    } else {
        for (size_t i = 0; i < vert.size(); ++i)
            if (vert[(start_vertex + vert.size() - i) % vert.size()] != vert_other[i])
                return false;
        return true;
    }
    return false;
}

bool Polygon::isConvex() const {
    bool flag_neg = false;
    bool flag_pos = false;
    for (size_t i = 0; i < verticesCount(); ++i) {
        size_t j = (i + 1) % verticesCount();
        size_t k = (i + 2) % verticesCount();
        double z = (vertices[j].x - vertices[i].x) * (vertices[k].y - vertices[j].y) 
                - (vertices[j].y - vertices[i].y) * (vertices[k].x - vertices[j].x);
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

void Polygon::rotate(const Point& _center, double angle) {
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

void Polygon::add_vertexes() {};

template<typename Head, typename... Tail>
void Polygon::add_vertexes(const Head& head, const Tail&... tail) {
    vertices.push_back(head);
    add_vertexes(tail...);
}

double Polygon::perimeter() const {
    double res = bandistance(vertices[0], vertices[vertices.size() - 1]);
    for (size_t i = 0; i < vertices.size() - 1; ++i)
        res += bandistance(vertices[i], vertices[i + 1]);
    return res;
}

double Polygon::area() const {
    double res = (vertices[0].x + vertices[vertices.size() - 1].x) * (vertices[0].y - vertices[vertices.size() - 1].y);
    for (size_t i = 0; i < vertices.size() - 1; ++i)
        res += (vertices[i].x + vertices[i + 1].x) * (vertices[i + 1].y - vertices[i].y);
    res /= 2;
    return abs(res);
}

bool Polygon::containsPoint(Point point) const {
    vector<Point> vert = getVertices();
    double sum = 0;
    for (size_t i = 0; i < vert.size(); ++i) {
        Point p = vert[i];
        Point q = vert[(i + 1) % vert.size()];
        Line first(p, point);
        Line second(q, point);
        double ang_cos = (first.a * second.a + first.b * second.b)
        / (sqrt(first.a * first.a + first.b * first.b) 
        * sqrt(second.a * second.a + second.b * second.b));
        double angle = acos(ang_cos);
        sum += angle;
    }
    return areSame(sum, 2. * M_PI);
};

class Triangle : public Polygon {
public:
    Triangle(Point a, Point b, Point c) : Polygon(a, b, c) {};
    Circle circumscribedCircle() const;
    Circle inscribedCircle() const;

    Point centroid() const;
    Point orthocenter() const;

    Line EulerLine() const;
    Circle ninePointsCircle() const;
};

Line Triangle::EulerLine() const {
    return Line(orthocenter(), circumscribedCircle().center());
}

Circle Triangle::ninePointsCircle() const {
    double centx =  (orthocenter().x + circumscribedCircle().center().x) / 2;
    double centy =  (orthocenter().y + circumscribedCircle().center().y) / 2;
    Point center(centx, centy);
    centx = (getVertices()[0].x + getVertices()[1].x) / 2;
    centy = (getVertices()[0].y + getVertices()[1].y) / 2;
    Point circ_point(centx, centy);
    double rad = bandistance(center, circ_point);
    return Circle(center, rad);
}

Point Triangle::centroid() const {
    vector<Point> verts = getVertices();
    Point mid1 = Point((verts[0].x + verts[1].x) / 2, (verts[0].y + verts[1].y) / 2);
    Point mid2 = Point((verts[1].x + verts[2].x) / 2, (verts[1].y + verts[2].y) / 2);
    Line median1 = Line(verts[2], mid1);
    Line median2 = Line(verts[0], mid2);
    return lineIntersection(median1, median2);
}

Point Triangle::orthocenter() const {
    vector<Point> vert = getVertices();

    double centx
    = (vert[2].x * vert[0].x + vert[1].y * vert[1].y) * vert[0].y
    + (vert[1].x * vert[2].x + vert[0].y * vert[0].y) * vert[2].y
    + (vert[0].x * vert[1].x + vert[2].y * vert[2].y) * vert[1].y
    - (vert[2].x * vert[0].x + vert[1].y * vert[1].y) * vert[2].y
    - (vert[1].x * vert[2].x + vert[0].y * vert[0].y) * vert[1].y
    - (vert[0].x * vert[1].x + vert[2].y * vert[2].y) * vert[0].y;

    double centy
    = (vert[0].x * vert[0].x + vert[1].y * vert[2].y) * vert[1].x
    + (vert[2].x * vert[2].x + vert[0].y * vert[1].y) * vert[0].x
    + (vert[1].x * vert[1].x + vert[2].y * vert[0].y) * vert[2].x
    - (vert[0].x * vert[0].x + vert[1].y * vert[2].y) * vert[2].x
    - (vert[2].x * vert[2].x + vert[0].y * vert[1].y) * vert[1].x
    - (vert[1].x * vert[1].x + vert[2].y * vert[0].y) * vert[0].x;

    double denom = vert[0].x * vert[1].y + vert[2].x * vert[0].y + vert[1].x * vert[2].y
    - vert[2].x * vert[1].y - vert[1].x * vert[0].y - vert[0].x * vert[2].y;

    centx /= denom;
    centy /= denom;

    return Point(centx, centy);
}

Circle Triangle::inscribedCircle() const {
    vector<Point> vert = getVertices();
    double a = bandistance(vert[1], vert[2]);
    double b = bandistance(vert[0], vert[2]);
    double c = bandistance(vert[0], vert[1]);
    double centx = (a * vert[0].x + b * vert[1].x + c * vert[2].x) / (a + b + c);
    double centy = (a * vert[0].y + b * vert[1].y + c * vert[2].y) / (a + b + c);
    Point center(centx, centy);
    double rad = 2 * area() / (a + b + c);
    return Circle(center, rad);
};

Circle Triangle::circumscribedCircle() const {
    vector<Point> vert = getVertices();
    double centx 
    = (vert[0].x * vert[0].x + vert[0].y * vert[0].y) * vert[1].y 
    + (vert[1].x * vert[1].x + vert[1].y * vert[1].y) * vert[2].y 
    + (vert[2].x * vert[2].x + vert[2].y * vert[2].y) * vert[0].y
    - (vert[0].x * vert[0].x + vert[0].y * vert[0].y) * vert[2].y 
    - (vert[1].x * vert[1].x + vert[1].y * vert[1].y) * vert[0].y 
    - (vert[2].x * vert[2].x + vert[2].y * vert[2].y) * vert[1].y;

    double centy 
    = (vert[0].x * vert[0].x + vert[0].y * vert[0].y) * vert[2].x
    + (vert[1].x * vert[1].x + vert[1].y * vert[1].y) * vert[0].x 
    + (vert[2].x * vert[2].x + vert[2].y * vert[2].y) * vert[1].x
    - (vert[0].x * vert[0].x + vert[0].y * vert[0].y) * vert[1].x 
    - (vert[1].x * vert[1].x + vert[1].y * vert[1].y) * vert[2].x 
    - (vert[2].x * vert[2].x + vert[2].y * vert[2].y) * vert[0].x;

    double denom = 2 * (vert[0].x * vert[1].y + vert[2].x * vert[0].y + vert[1].x * vert[2].y
    - vert[2].x * vert[1].y - vert[1].x * vert[0].y - vert[0].x * vert[2].y);

    centx /= denom;
    centy /= denom;

    Point center(centx, centy);
    double rad = bandistance(center, vert[0]);

    return Circle(center, rad);
}

class Rectangle : public Polygon {
public:
    Rectangle(Point x, Point y, double coef);
    Rectangle(Point a, Point b, Point c, Point d) : Polygon(a, b, c, d) {};
    Point center() const;
    pair<Line, Line> diagonals() const;

private:

};

Rectangle::Rectangle(Point x, Point y, double coef) {
    Point cent((x.x + y.x) / 2, (x.y + y.y) / 2);
    coef = max(coef, 1.0 / coef);
    double diag_len = bandistance(x, y);
    double width = diag_len / sqrt(1.0 + coef * coef);
    
    Line diag(x, y);
    Point p(diag.a, diag.b);
    Point diag_vec(y.x - x.x, y.y - x.y);
    Point height_vec(diag_vec.x / (1.0 + coef * coef), diag_vec.y / (1.0 + coef * coef));
    Point h(x.x + height_vec.x, x.y + height_vec.y);

    p = Point(p.x / sqrt(p.x * p.x + p.y * p.y), p.y / sqrt(p.x * p.x + p.y * p.y));
    double case_true_x = h.x + p.x * (diag_len * coef / (coef * coef + 1.0));
    double case_true_y = h.y + p.y * (diag_len * coef / (coef * coef + 1.0));

    double case_false_x = h.x - p.x * (diag_len * coef / (coef * coef + 1.0));
    double case_false_y = h.y - p.y * (diag_len * coef / (coef * coef + 1.0));

    Point case_true(case_true_x, case_true_y);
    Point case_false(case_false_x, case_false_y);

    Point a;
    if (areSame(bandistance(case_true, x), width) && !areSame(p.x * diag_vec.y - p.y * diag_vec.x, 0.))
        a = case_true;
    Point b(a.x + (cent.x - a.x) * 2, a.y + (cent.y - a.y) * 2);
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
    Square(Point a, Point b) : Rectangle(a, b, 1) {};
    Circle circumscribedCircle() const;
    Circle inscribedCircle() const;
};

Circle Square::inscribedCircle() const {
    return Circle(center(), bandistance(getVertices()[0], getVertices()[1]) / 2);
}

Circle Square::circumscribedCircle() const {
    return Circle(center(), bandistance(center(), getVertices()[0]));
}

int main() {
    Point a(0, 0);
    Point b(2, 1);
    Point c(6, 3);
    Rectangle r1(a, b, 2);
    Rectangle r2(b, c, 2);
    cout << (r1.isSimilarTo(r2) ? 1 : 0) << endl;
    system("pause");
    return 0;
}