# Using Graphviz

```bash
echo 'graph { 0--1; 0--2; 0--11; 1--3; 1--8; 2--3; 2--9; 11--8; 11--9; 8--5; 9--5; 3--5 }' | dot -Tsvg > output-dot1.svg
echo 'graph { 0--1; 0--2; 0--11; 1--3; 1--8; 2--3; 2--9; 11--8; 11--9; 5--8; 9--5; 5--3 }' | dot -Tsvg > output-dot2.svg
echo 'graph { 0--1; 0--2; 0--11; 1--3; 1--8; 2--3; 2--9; 11--8; 11--9; 5--8; 9--5; 5--3 }' | neato -Tsvg > output-neato.svg
echo 'graph { 0--1; 0--2; 0--11; 1--3; 1--8; 2--3; 2--9; 11--8; 11--9; 5--8; 9--5; 5--3 }' | sfdp -Tsvg > output-sfdp.svg
echo 'graph { 0--1; 0--2; 0--11; 1--3; 1--8; 2--3; 2--9; 11--8; 11--9; 5--8; 9--5; 5--3 }' | circo -Tsvg > output-circo.svg
```

