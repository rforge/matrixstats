data(loexp.example)

plot(ex1, pch=19, type="b")

# 50th expectile (median)
lines(loexp(ex1, alpha=0.5)$outy, lwd=4, col="blue")

# 0.5th expectile (baseline)
lines(loexp(ex1, alpha=0.005)$outy, lwd=4, col="red")

plot(ex2.y, pch=19, type="b")
# give weight=0 for 0s
lines(loexp(ex2.y, w=ex2.w, alpha=0.005)$outy, lwd=4, col="blue")
