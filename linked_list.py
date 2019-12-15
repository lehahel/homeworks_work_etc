class Node:
    def __init__(self, _key, _next):
        self.key = _key
        self.next = _next

    def __del__(self):
        del self.next

    def insert(self, _value):
        cur = self
        while cur.next is not None:
            cur = cur.next
        cur.next = _value

    def erase(self, _key):
        cur = self
        while cur is not None:
            if cur.key == _key:
                old_data = cur
                cur = old_data.next
                old_data.next = None
                del old_data
                return True
            cur = cur.next
        return False

    def print(self, node):
        if node is not None:
            print(node.key, end=' ')
            node.print(node.next)

    def reverse_print(self, node):
        if node is not None:
            node.reverse_print(node.next)
            print(node.key, end=' ')

    def perform_bubble_sort_step(self, node):
        if node.next is not None:
            if node.key.upper() > node.next.key.upper():
                node.key, node.next.key = node.next.key, node.key
            node.perform_bubble_sort_step(node.next)


class SingleLinkedList:
    def __init__(self):
        self.head = None
        self.size = 0

    def __del__(self):
        del self.head

    def insert(self, _key):
        if not self.head:
            self.head = Node(_key, None)
            self.size += 1
        else:
            Node.insert(self.head, Node(_key, None))
            self.size += 1

    def erase(self, _key):
        if self.size == 0:
            return False
        if self.head.key == _key:
            self.head = self.head.next
            self.size -= 1
        elif Node.erase(self.head, _key):
            self.size -= 1
            return True
        return False

    def print(self):
        Node.print(self.head, self.head)
        print()

    def reverse_print(self):
        Node.reverse_print(self.head, self.head)
        print()

    def sort(self):
        if self.size > 1:
            for i in range(self.size):
                Node.perform_bubble_sort_step(self.head, self.head)


string = input()
cur_word = ""
arr = SingleLinkedList()
for i in range(len(string)):
    if 'a' <= string[i] <= 'z' or 'A' <= string[i] <= 'Z':
        cur_word += string[i]
    else:
        arr.insert(cur_word)
        cur_word = ""

if cur_word != "":
    arr.insert(cur_word)

arr.sort()
arr.print()
arr.reverse_print()
